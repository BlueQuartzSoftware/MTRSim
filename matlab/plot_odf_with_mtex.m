function plot_odf_with_mtex()
    % plot_odf_with_mtex
    %
    % Build an MTEX ODF from the masked phase-1 voxels of
    % data/real_world_microtexture_data.dream3d (NOT from our binned
    % calc_ODF.m output) and produce three figures:
    %
    %   1. Standard HCP pole-figure montage: {0001}, {10-10}, {11-20}
    %   2. Phi2 sections via plot(odf)
    %   3. 3-D Euler-space rendering via plot3d(odf)
    %
    % This is the natural MTEX path: feed it raw EBSD orientations and
    % let it kernel-density-estimate the ODF on its own. Result is a
    % texture visualisation independent of our calc_ODF.m bin pipeline,
    % suitable for visual sanity-check and report figures.
    %
    % Run:
    %   matlab -batch "addpath('/.../MTRSim/matlab'); plot_odf_with_mtex()"

    % --- Self-locate repo and ensure MTEX is on the path ---
    repoRoot = fileparts(fileparts(mfilename('fullpath')));

    % MTEX must be initialised. Try the standard MTEX startup if it isn't
    % already on the path. Search a few well-known install locations.
    if exist('crystalSymmetry', 'file') ~= 2
        if exist('startup_mtex', 'file') == 2
            startup_mtex;
        else
            mtexCandidates = { ...
                fullfile(repoRoot, '..', 'mtex-6.1.0'), ...
                fullfile(repoRoot, '..', '..', 'mtex-6.1.0'), ...
                '/Users/mjackson/Workspace7/mtex-6.1.0' ...
            };
            mtexFound = false;
            for k = 1:numel(mtexCandidates)
                startupPath = fullfile(mtexCandidates{k}, 'startup_mtex.m');
                if exist(startupPath, 'file') == 2
                    fprintf('Initialising MTEX from %s\n', mtexCandidates{k});
                    addpath(mtexCandidates{k});
                    run(startupPath);
                    mtexFound = true;
                    break;
                end
            end
            if ~mtexFound
                error(['MTEX is not on the MATLAB path and no startup_mtex.m was ' ...
                       'found in the usual locations. Run startup_mtex.m from your ' ...
                       'MTEX install manually before calling this function.']);
            end
        end
    end

    outDir = fullfile(repoRoot, 'output', 'mtex_figures');
    if ~exist(outDir, 'dir'); mkdir(outDir); end

    % --- Load the EBSD .dream3d and filter to contributing voxels ---
    inputFile = fullfile(repoRoot, 'data', 'real_world_microtexture_data.dream3d');
    if ~exist(inputFile, 'file')
        error('Input file not found: %s', inputFile);
    end

    fprintf('Reading EBSD data from %s...\n', inputFile);
    eulers  = h5read(inputFile, '/DataStructure/16_Micro_Res/CellData/EulerAngles');
    phases  = h5read(inputFile, '/DataStructure/16_Micro_Res/CellData/Phases');
    mask    = h5read(inputFile, '/DataStructure/16_Micro_Res/CellData/Mask');

    eulers = reshape(eulers, 3, []);
    phases = phases(:);
    mask   = mask(:);
    valid  = (phases == 1) & (mask == 1);

    phi1 = double(eulers(1, valid))';
    PHI  = double(eulers(2, valid))';
    phi2 = double(eulers(3, valid))';

    fprintf('  Contributing voxels: %d\n', numel(phi1));

    % --- Build MTEX orientation set ---
    % Hexagonal 6/mmm with approximate Ti alpha lattice constants.
    cs = crystalSymmetry('6/mmm', [3 3 4.7], 'mineral', 'Titanium');
    ori = orientation.byEuler(phi1, PHI, phi2, cs);

    fprintf('Building MTEX ODF from %d orientations (kernel halfwidth = 5 deg)...\n', numel(ori));
    odf = calcDensity(ori, 'halfwidth', 5*degree);
    fprintf('  Texture index (squared norm): %.3f\n', norm(odf)^2);
    fprintf('  Max ODF value (MUD): %.3f\n', max(odf));

    % --- Figure 1: HCP pole-figure montage ---
    fprintf('\n[1/3] Pole figures: {0001}, {10-10}, {11-20}\n');
    h = [Miller(0, 0, 0, 1, cs), Miller(1, 0, -1, 0, cs), Miller(1, 1, -2, 0, cs)];
    figure('Position', [100 100 1200 400], 'Name', 'HCP Pole Figures (MTEX from raw EBSD)');
    plotPDF(odf, h, 'antipodal');
    polePath = fullfile(outDir, 'pole_figures_0001_1010_1120.png');
    saveas(gcf, polePath);
    fprintf('  saved %s\n', polePath);

    % --- Figure 2: ODF phi2 sections (the "textbook" ODF view) ---
    fprintf('\n[2/3] ODF phi2 sections via plot(odf)\n');
    figure('Position', [100 100 1200 800], 'Name', 'ODF phi2 Sections (MTEX)');
    plot(odf);
    sectPath = fullfile(outDir, 'odf_phi2_sections.png');
    saveas(gcf, sectPath);
    fprintf('  saved %s\n', sectPath);

    % --- Figure 3: 3D Euler-space rendering ---
    fprintf('\n[3/3] 3D ODF rendering via plot3d(odf)\n');
    figure('Position', [100 100 900 800], 'Name', 'ODF 3D Euler-space (MTEX)');
    plot3d(odf);
    threeDPath = fullfile(outDir, 'odf_3d_euler_space.png');
    saveas(gcf, threeDPath);
    fprintf('  saved %s\n', threeDPath);

    fprintf('\nDone. Figures saved to %s\n', outDir);
end
