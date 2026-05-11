function run_validation()
    % run_validation
    %
    % Top-level dispatcher for the MATLAB-side ODF validation reference
    % generation. Run from anywhere as:
    %   matlab -batch "addpath('/abs/path/to/MTRSim/matlab'); run_validation()"
    % or from the repo root as:
    %   matlab -batch "addpath('matlab'); run_validation()"
    %
    % Phase 1 (this file): produces data/calc_odf_reference_targeted.h5 from
    % a hardcoded 12-orientation HCP fixture for bin-by-bin diff against the
    % C++ ComputeODFFilter output.
    %
    % Phase 2 (run_realistic_validation, not yet implemented) will validate
    % against the masked phase-1 voxels in
    % data/real_world_microtexture_data.dream3d.

    % Self-locate the repo root and ensure both matlab/ and data/ are on the
    % MATLAB path. calc_ODF.m loads blank_ODF.mat from data/ via `load`, so
    % data/ must be visible.
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(repoRoot, 'matlab'));
    addpath(fullfile(repoRoot, 'data'));

    run_targeted_validation();
    run_realistic_validation(repoRoot);
end

function run_realistic_validation(repoRoot)
    % Phase 2: validate against the masked phase-1 voxels of
    % data/real_world_microtexture_data.dream3d. The .dream3d file is
    % standard HDF5 with the DREAM3D-NX layout:
    %   /DataStructure/16_Micro_Res/CellData/EulerAngles  (1, 640, 640, 3) float32, radians
    %   /DataStructure/16_Micro_Res/CellData/Phases       (1, 640, 640, 1) int32  (1 = Titanium Alpha)
    %   /DataStructure/16_Micro_Res/CellData/Mask         (1, 640, 640, 1) uint8  (1 = valid)
    %   /DataStructure/16_Micro_Res/CellEnsembleData/CrystalStructures (2,1) uint32 (phase 1 = 0 = Hexagonal_High)
    %
    % We extract only the voxels where mask==1 AND phase==1 (~357k of 409.6k).

    inputFile = fullfile(repoRoot, 'data', 'real_world_microtexture_data.dream3d');
    if ~exist(inputFile, 'file')
        fprintf('Skipping run_realistic_validation: %s not found\n', inputFile);
        return;
    end

    fprintf('Reading EBSD data from %s...\n', inputFile);
    eulers  = h5read(inputFile, '/DataStructure/16_Micro_Res/CellData/EulerAngles'); % (3, 640, 640, 1) float32 (h5read transposes)
    phases  = h5read(inputFile, '/DataStructure/16_Micro_Res/CellData/Phases');      % (1, 640, 640, 1) int32
    mask    = h5read(inputFile, '/DataStructure/16_Micro_Res/CellData/Mask');        % (1, 640, 640, 1) uint8

    % h5read transposes dimensions vs the on-disk shape. EulerAngles becomes
    % (3, 640, 640, 1) — flatten across spatial axes.
    eulers  = reshape(eulers, 3, []);     % (3, 409600)
    phases  = phases(:);                  % (409600, 1)
    mask    = mask(:);                    % (409600, 1)

    % Build the contributing-voxel mask: phase==1 (HCP) AND mask==1.
    valid = (phases == 1) & (mask == 1);
    fprintf('  Total voxels:       %d\n', numel(phases));
    fprintf('  Phase==1 voxels:    %d\n', sum(phases == 1));
    fprintf('  Mask==1 voxels:     %d\n', sum(mask == 1));
    fprintf('  Contributing (both): %d\n', sum(valid));

    phi1 = double(eulers(1, valid))';    % column vector, double precision
    PHI  = double(eulers(2, valid))';
    phi2 = double(eulers(3, valid))';

    fprintf('Calling calc_ODF on %d orientations (this will take a few minutes)...\n', numel(phi1));
    [ODFval, ~, phi1_bins, PHI_bins, phi2_bins] = calc_ODF(phi1, PHI, phi2);

    out = fullfile(repoRoot, 'data', 'calc_odf_reference_realistic.h5');
    if exist(out, 'file'); delete(out); end

    h5create(out, '/ODF_best/num_components', [1], 'Datatype', 'int64');
    h5write(out,  '/ODF_best/num_components', int64(1));

    h5create(out, '/ODF_best/component_0/ODFval', size(ODFval), 'Datatype', 'double');
    h5write(out,  '/ODF_best/component_0/ODFval', ODFval);

    h5create(out, '/ODF_best/component_0/phi1_bins', size(phi1_bins), 'Datatype', 'double');
    h5write(out,  '/ODF_best/component_0/phi1_bins', phi1_bins);

    h5create(out, '/ODF_best/component_0/PHI_bins', size(PHI_bins), 'Datatype', 'double');
    h5write(out,  '/ODF_best/component_0/PHI_bins', PHI_bins);

    h5create(out, '/ODF_best/component_0/phi2_bins', size(phi2_bins), 'Datatype', 'double');
    h5write(out,  '/ODF_best/component_0/phi2_bins', phi2_bins);

    fprintf('\nWrote realistic reference ODF to %s\n', out);
    fprintf('  ODFval size: %d, sum: %.6f, max: %.6e, nonzero: %d\n', ...
            numel(ODFval), sum(ODFval(:)), max(ODFval(:)), nnz(ODFval));
end

function run_targeted_validation()
    % Hardcoded 12-orientation HCP fixture.
    %
    % FP-precision-safe fixture design (important):
    %   Every input value is of the form (5 deg * N + 2.5 deg) -- i.e. it sits
    %   at the EXACT MID-POINT of its 5 deg bin. Hex 6/mmm symmetry operators
    %   are at multiples of 60 deg (z-rot) and 30 deg (secondary 2-fold), so
    %   every symmetric variant of a mid-bin input is also mid-bin in all
    %   three Bunge axes. That gives 0.5 step-units (= 2.5 deg) of margin from
    %   the nearest bin boundary on every axis -- many orders of magnitude
    %   beyond the ulp-level FP precision of atan2/floor/fix at exact-boundary
    %   inputs. As a result, MATLAB / Python / C++ all agree bin-for-bin.
    %
    %   See test/ComputeODFTest.cpp for the full discussion of why the earlier
    %   exact-boundary fixture (phi1=45, PHI=0/90/180, phi1/phi2=359, etc.)
    %   produced ~1851 bins of disagreement at ulp-level FP precision pathology.
    %
    % This table is byte-identical to the comment + array in
    % test/ComputeODFTest.cpp. A future reader should be able to diff the two
    % and see they're the same.
    %
    %  # | (phi1 deg, PHI deg, phi2 deg) | bin (5 deg spacing)   | Why
    % ---+--------------------------------+-----------------------+----------------------------------------------
    %   1| ( 12.5,    12.5,   12.5)       | (2,  2,  2)           | Pure mid-bin; baseline
    %   2| ( 47.5,    27.5,   92.5)       | (9,  5,  18)          | Generic interior, no special structure
    %   3| (137.5,    67.5,  217.5)       | (27, 13, 43)          | Multi-decimal-bin coverage; phi1/phi2 > 90 deg
    %   4| (  2.5,     2.5,    2.5)       | (0,  0,  0)           | All near zero; tests near-PHI=0 handling
    %   5| ( 47.5,     2.5,   32.5)       | (9,  0,  6)           | PHI very small with non-trivial phi1/phi2
    %   6| ( 62.5,     2.5,   92.5)       | (12, 0,  18)          | PHI very small, generic phi1/phi2
    %   7| ( 47.5,    92.5,   32.5)       | (9,  18, 6)           | PHI just above pi/2 (equatorial regime)
    %   8| ( 47.5,   177.5,   32.5)       | (9,  35, 6)           | PHI near pi (upper-pole regime)
    %   9| ( 47.5,    87.5,   32.5)       | (9,  17, 6)           | PHI just below pi/2 (complement of #7)
    %  10| (357.5,    32.5,   32.5)       | (71, 6,  6)           | phi1 near 360 deg (wrap regime)
    %  11| ( 47.5,    32.5,  357.5)       | (9,  6,  71)          | phi2 near 360 deg (wrap regime)
    %  12| (357.5,    87.5,  357.5)       | (71, 17, 71)          | All three near upper boundaries simultaneously
    eulersDeg = [
         12.5,  12.5,   12.5;
         47.5,  27.5,   92.5;
        137.5,  67.5,  217.5;
          2.5,   2.5,    2.5;
         47.5,   2.5,   32.5;
         62.5,   2.5,   92.5;
         47.5,  92.5,   32.5;
         47.5, 177.5,   32.5;
         47.5,  87.5,   32.5;
        357.5,  32.5,   32.5;
         47.5,  32.5,  357.5;
        357.5,  87.5,  357.5
    ];

    eulersRad = deg2rad(eulersDeg);
    phi1 = eulersRad(:, 1);
    PHI  = eulersRad(:, 2);
    phi2 = eulersRad(:, 3);

    % calc_ODF.m signature: [ODFval, ODFbins, phi1_bins, PHI_bins, phi2_bins] = calc_ODF(phi1, PHI, phi2)
    % It hardcodes degree_spacing=5 (=> 72 x 36 x 72 = 186624 bins),
    % applies HCP symmetry expansion via symmetric_euler_angles.m (12 ops),
    % and applies the tri-linear smoothing kernel (smooth_ODF=true) by default.
    % Our C++ ComputeODFFilter defaults apply_smoothing=true and bin_size_deg=5
    % to match.
    fprintf('Calling calc_ODF on %d orientations...\n', size(eulersDeg, 1));
    [ODFval, ~, phi1_bins, PHI_bins, phi2_bins] = calc_ODF(phi1, PHI, phi2);

    % Reshape bin-edge vectors to row vectors so h5write writes them as 1-D
    % datasets of the right length (73, 37, 73) without dimension-mismatch
    % surprises across MATLAB versions.
    phi1_bins = phi1_bins(:)';
    PHI_bins  = PHI_bins(:)';
    phi2_bins = phi2_bins(:)';

    % Ensure ODFval is a flat column vector of length 186624 so h5write
    % stores it as a 1-D dataset that mtrsim::readODFComponents reads as a
    % flat std::vector<double>.
    ODFval = ODFval(:);

    % Output path: matlab/run_validation.m -> matlab/.. -> data/calc_odf_reference_targeted.h5
    out = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'calc_odf_reference_targeted.h5');
    if exist(out, 'file')
        delete(out);
    end

    % Layout matches src/LibMTRSim/ODFFileIO.{hpp,cpp}'s readODFMetadata /
    % readODFComponents expectation:
    %   /ODF_best/num_components            (int64 scalar)
    %   /ODF_best/component_0/ODFval        (float64, length 186624)
    %   /ODF_best/component_0/phi1_bins     (float64, edges in radians, length 73)
    %   /ODF_best/component_0/PHI_bins      (float64, edges in radians, length 37)
    %   /ODF_best/component_0/phi2_bins     (float64, edges in radians, length 73)
    %
    % Plus a fixture-version field so the C++ test can detect a stale reference
    % HDF5 (mismatched orientations between this MATLAB script and the test) and
    % skip with a clear regeneration message instead of producing a confusing
    % bin-by-bin diff:
    %   /ODF_best/fixture_version           (int64 scalar; bump on fixture change)
    %
    % Bump this value (and the matching k_TargetedFixtureVersion in
    % test/ComputeODFTest.cpp) whenever the orientation table above changes.
    %   v1: original boundary-stress fixture (45 deg, 90 deg, 180 deg, 359 deg, ...)
    %   v2: FP-precision-safe mid-bin fixture (every input is N*5 + 2.5 deg)
    fixture_version = int64(2);
    h5create(out, '/ODF_best/fixture_version', 1, 'Datatype', 'int64');
    h5write(out,  '/ODF_best/fixture_version', fixture_version);

    h5create(out, '/ODF_best/num_components', 1, 'Datatype', 'int64');
    h5write(out,  '/ODF_best/num_components', int64(1));

    h5create(out, '/ODF_best/component_0/ODFval', size(ODFval), 'Datatype', 'double');
    h5write(out,  '/ODF_best/component_0/ODFval', ODFval);

    h5create(out, '/ODF_best/component_0/phi1_bins', size(phi1_bins), 'Datatype', 'double');
    h5write(out,  '/ODF_best/component_0/phi1_bins', phi1_bins);

    h5create(out, '/ODF_best/component_0/PHI_bins', size(PHI_bins), 'Datatype', 'double');
    h5write(out,  '/ODF_best/component_0/PHI_bins', PHI_bins);

    h5create(out, '/ODF_best/component_0/phi2_bins', size(phi2_bins), 'Datatype', 'double');
    h5write(out,  '/ODF_best/component_0/phi2_bins', phi2_bins);

    fprintf('Wrote targeted reference ODF to %s\n', out);
    fprintf('  ODFval size: %d, sum: %.6f, max: %.6f\n', numel(ODFval), sum(ODFval(:)), max(ODFval(:)));
end
