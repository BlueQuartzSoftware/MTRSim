function run_simulation_from_json(config_path, output_dir)
% run_simulation_from_json  Run MTR simulation driven by a JSON config file.
%
%   run_simulation_from_json(CONFIG_PATH, OUTPUT_DIR)
%
%   This wrapper reads the same JSON configuration files used by the C++
%   mtrsim executable and produces output in an identical CSV format so that
%   the two implementations can be compared.
%
%   Outputs (written to OUTPUT_DIR):
%       matlab_sim_results.csv   – x,y,z,phi1,PHI,phi2,mtr_index
%       matlab_sim_IPF_map.jpg   – Inverse Pole Figure colour map
%       matlab_sim_PGRF.csv      – x,y,z,mtr_index  (PGRF-only results)
%
%   Example (command-line):
%       matlab -nodisplay -nosplash -batch \
%           "run_simulation_from_json('configs/default.json','output/matlab')"

%% ---- resolve paths ---------------------------------------------------------
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);          % one level up from matlab/

if nargin < 2
    output_dir = fullfile(project_root, 'output', 'matlab');
end

if ~isfolder(output_dir)
    mkdir(output_dir);
end

%% ---- read JSON config ------------------------------------------------------
config_text = fileread(config_path);
cfg = jsondecode(config_text);

% Map JSON keys to local variables (with defaults matching SimulationParams.hpp)
xlen = getfield_or(cfg, 'xLen', 1.5 * 25.4);
ylen = getfield_or(cfg, 'yLen', 0.5 * 25.4);
zlen = getfield_or(cfg, 'zLen', 0.0);
dx   = getfield_or(cfg, 'dx',  0.02);
dy   = getfield_or(cfg, 'dy',  0.02);
dz   = getfield_or(cfg, 'dz',  0.02);

volume_fractions = getfield_or(cfg, 'volumeFractions', [0.30, 0.35, 0.35]);
volume_fractions = volume_fractions(:)';  % ensure row vector

theta_list = getfield_or(cfg, 'thetaList', [0.10, 0.45, 0.10; 0.08, 0.37, 0.08]);
if iscell(theta_list)
    theta_list = cell2mat(theta_list);
end

nugvar = getfield_or(cfg, 'nuggetVariance', [0.67, 0.71, 0.72]);
nugvar = nugvar(:)';

seed_val = getfield_or(cfg, 'seed', 0);

%% ---- seed the RNG ----------------------------------------------------------
if seed_val == 0
    rng('shuffle');
    fprintf('[MATLAB] Using random seed (shuffle)\n');
else
    rng(seed_val);
    fprintf('[MATLAB] Using fixed seed: %d\n', seed_val);
end

%% ---- compute grid -----------------------------------------------------------
nx = round(xlen / dx);
ny = round(ylen / dy);
nz = max(round(zlen / dz), 1);
N  = nx * ny * nz;
P  = volume_fractions;

fprintf('[MATLAB] Grid: nx=%d  ny=%d  nz=%d  N=%d\n', nx, ny, nz, N);
fprintf('[MATLAB]   xLen=%.3f mm  yLen=%.3f mm  zLen=%.3f mm\n', xlen, ylen, zlen);
fprintf('[MATLAB]   dx=%.4f  dy=%.4f  dz=%.4f [mm]\n', dx, dy, dz);

%% ---- build spatial coordinates (same order as C++) --------------------------
% z outer -> x middle -> y inner
s_sim = zeros(N, 3);
k = 0;
for zix = 1:nz
    z_tmp = zix * dz;
    for j = 1:nx
        x_tmp = j * dx;
        for i = 1:ny
            y_tmp = i * dy;
            k = k + 1;
            s_sim(k, :) = [x_tmp, y_tmp, z_tmp];
        end
    end
end

%% ---- run PGRF simulation ---------------------------------------------------
fprintf('[MATLAB] Running PGRF simulation...\n');
tic;
[MTR_IX_sim, Y_sim] = PGRF_simulation(P, theta_list, nugvar, dx, dy, dz, nx, ny, nz);
pgrf_time = toc;
fprintf('[MATLAB] PGRF complete (%.2f s)\n', pgrf_time);

% Report actual volume fractions
num_components = numel(P);
fprintf('[MATLAB] Volume fractions (target -> actual):\n');
for j = 1:num_components
    actual = sum(MTR_IX_sim == j) / N;
    fprintf('  Component %d: %.3f -> %.3f\n', j, P(j), actual);
end

%% ---- write PGRF-only CSV ---------------------------------------------------
pgrf_csv_path = fullfile(output_dir, 'matlab_sim_PGRF.csv');
fprintf('[MATLAB] Writing PGRF CSV: %s\n', pgrf_csv_path);
fid = fopen(pgrf_csv_path, 'w');
fprintf(fid, 'x,y,z,mtr_index\n');
for i = 1:N
    fprintf(fid, '%.6f,%.6f,%.6f,%d\n', s_sim(i,1), s_sim(i,2), s_sim(i,3), MTR_IX_sim(i));
end
fclose(fid);

%% ---- load ODF data ---------------------------------------------------------
% The .mat files are in the data/ directory relative to project root
data_dir = fullfile(project_root, 'data');
fprintf('[MATLAB] Loading ODF data from: %s\n', data_dir);

uniform_file = fullfile(data_dir, 'uniformODF.mat');
odf_file = fullfile(data_dir, 'simulation_ODF.mat');

if ~isfile(uniform_file) || ~isfile(odf_file)
    warning('ODF .mat files not found in %s. Skipping orientation sampling.', data_dir);
    % Write PGRF-only results
    csv_path = fullfile(output_dir, 'matlab_sim_results.csv');
    fid = fopen(csv_path, 'w');
    fprintf(fid, 'x,y,z,phi1,PHI,phi2,mtr_index\n');
    for i = 1:N
        fprintf(fid, '%.6f,%.6f,%.6f,0.000000,0.000000,0.000000,%d\n', ...
            s_sim(i,1), s_sim(i,2), s_sim(i,3), MTR_IX_sim(i));
    end
    fclose(fid);
    fprintf('[MATLAB] Done (PGRF only, no ODF data).\n');
    return;
end

load(uniform_file, 'uniformODF');
tmp = load(odf_file, 'ODF_best');
ODF_best = tmp.ODF_best;

% Normalise ODF values to integrate to 1 (matches C++ preprocessing)
for j = 1:numel(ODF_best)
    ODF_best(j).ODFval = ODF_best(j).ODFval / sum(ODF_best(j).ODFval);
end

%% ---- sample orientations ---------------------------------------------------
fprintf('[MATLAB] Sampling orientations (N=%d)...\n', N);
tic;

phi1_vec_sim = zeros(N, 1);
PHI_vec_sim  = zeros(N, 1);
phi2_vec_sim = zeros(N, 1);

% Pre-sample N orientations per component (batched, matches C++ approach)
phi1_a = zeros(N, numel(ODF_best));
PHI_a  = zeros(N, numel(ODF_best));
phi2_a = zeros(N, numel(ODF_best));

for j = 1:numel(ODF_best)
    fprintf('  Component %d / %d\n', j, numel(ODF_best));
    [phi1_accept, PHI_accept, phi2_accept] = sample_N_orientations_from_ODF(N, ODF_best(j), uniformODF);
    phi1_a(:, j) = phi1_accept;
    PHI_a(:, j)  = PHI_accept;
    phi2_a(:, j) = phi2_accept;
end

% Assign each voxel its orientation based on MTR component
for i = 1:N
    phi1_vec_sim(i, 1) = phi1_a(i, MTR_IX_sim(i, 1));
    PHI_vec_sim(i, 1)  = PHI_a(i, MTR_IX_sim(i, 1));
    phi2_vec_sim(i, 1) = phi2_a(i, MTR_IX_sim(i, 1));
end

orient_time = toc;
fprintf('[MATLAB] Orientation sampling complete (%.2f s)\n', orient_time);

%% ---- write full results CSV ------------------------------------------------
csv_path = fullfile(output_dir, 'matlab_sim_results.csv');
fprintf('[MATLAB] Writing results CSV: %s\n', csv_path);
fid = fopen(csv_path, 'w');
fprintf(fid, 'x,y,z,phi1,PHI,phi2,mtr_index\n');
for i = 1:N
    fprintf(fid, '%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d\n', ...
        s_sim(i,1), s_sim(i,2), s_sim(i,3), ...
        phi1_vec_sim(i), PHI_vec_sim(i), phi2_vec_sim(i), ...
        MTR_IX_sim(i));
end
fclose(fid);

%% ---- write IPF map ----------------------------------------------------------
fprintf('[MATLAB] Generating IPF map...\n');
try
    a = view_IPF_map(s_sim, phi1_vec_sim, PHI_vec_sim, phi2_vec_sim);
    ipf_path = fullfile(output_dir, 'matlab_sim_IPF_map.jpg');
    imwrite(a, ipf_path, 'jpg');
    fprintf('[MATLAB] IPF map written: %s\n', ipf_path);
catch me
    fprintf('[MATLAB] IPF map generation failed: %s\n', me.message);
end

%% ---- done -------------------------------------------------------------------
total_time = pgrf_time + orient_time;
fprintf('[MATLAB] Total simulation time: %.2f s\n', total_time);
fprintf('[MATLAB] Done.\n');

end

%% ---- helper: safe field access with default ---------------------------------
function val = getfield_or(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end
