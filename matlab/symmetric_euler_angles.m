% symmetric_euler_angles.m

% Daniel M. Sparkman
% 02/14/2013
% Research

% symmetric_euler_angles.m

% function [X,Y] = generate_inverse_pole_figure(phi1,PHI,phi2)
function [phi1_f, PHI_f, phi2_f] = symmetric_euler_angles(phi1,PHI,phi2)

% hcp     fcc
crystal_system = 'hcp';


%%

% PHI = mod(PHI,pi/180*90);
% phi2 = mod(phi2,pi/180*60);

num_points = numel(phi1);

%% develop symmetry operations

all_symmetry = true;

if all_symmetry

    switch crystal_system
        case 'hcp'
            % number of symmetry operations
            N = 12;
            
            symmetry_operators_euler_angles = pi/180* [ 0 0  0;
                0 0 60;
                0 0 120;
                0 0 180;
                0 0 240;
                0 0 300;
                0 180 0;
                0 180 60;
                0 180 120;
                0 180 180;
                0 180 240;
                0 180 300];
            
        case 'fcc'
            % cubic
            % checked with MCS for max theta_mis = 62.37
            % number of symmetry operations
            N = 24;
            
            symmetry_operators_euler_angles = pi/180* [ 0 0  0;
                0 0 90;
                0 0 180;
                0 0 270;
                0 90 0;
                0 90 90;
                0 90 180;
                0 90 270;
                0 180 0;
                0 180 90;
                0 180 180;
                0 180 270;
                0 270 0;
                0 270 90;
                0 270 180;
                0 270 270;
                90 90 0;
                90 90 90;
                90 90 180;
                90 90 270;
                90 270 0;
                90 270 90;
                90 270 180;
                90 270 270];
    end

    phi1s = symmetry_operators_euler_angles(:,1)';
    PHIs = symmetry_operators_euler_angles(:,2)';
    phi2s = symmetry_operators_euler_angles(:,3)';
    
    % just the symmetry operators for one orientation
    phi1_symm = zeros(N,1);
    PHI_symm  = zeros(N,1);
    phi2_symm = zeros(N,1);
    
    % all symmetric orientations
    phi1_f = zeros(N*num_points,1);
    PHI_f  = zeros(N*num_points,1);
    phi2_f = zeros(N*num_points,1);
end

%% step through each orientation

c1 = cos(phi1);
c2 = cos(phi2);
C = cos(PHI);
s1 = sin(phi1);
s2 = sin(phi2);
S = sin(PHI);

% build rotation matrix for input orientations
g11 = c1.*c2-s1.*s2.*C;
g12 = s1.*c2+c1.*s2.*C;
g13 = s2.*S;
g21 = -c1.*s2-s1.*c2.*C;
g22 = -s1.*s2+c1.*c2.*C;
g23 = c2.*S;
g31 = s1.*S;
g32 = -c1.*S;
g33 = C;

% setup symmetry operators (each should be 1xN)
c1s = cos(phi1s);
c2s = cos(phi2s);
Cs  = cos(PHIs);
s1s = sin(phi1s);
s2s = sin(phi2s);
Ss  = sin(PHIs);

% build crystal symmetry operator rotation matrix (each should be 1xN)
O_crystal_11 = c1s.*c2s-s1s.*s2s.*Cs;
O_crystal_12 = s1s.*c2s+c1s.*s2s.*Cs;
O_crystal_13 = s2s.*Ss;
O_crystal_21 = -c1s.*s2s-s1s.*c2s.*Cs;
O_crystal_22 = -s1s.*s2s+c1s.*c2s.*Cs;
O_crystal_23 = c2s.*Ss;
O_crystal_31 = s1s.*Ss;
O_crystal_32 = -c1s.*Ss;
O_crystal_33 = Cs;


% build symmetry rotation matrices (each should be num_pointsxN)
% g_rot_symm_11 = g11 * O_crystal_11 + g21 * O_crystal_12 + g31 * O_crystal_13;
% g_rot_symm_12 = g12 * O_crystal_11 + g22 * O_crystal_12 + g32 * O_crystal_13;
g_rot_symm_13 = g13 * O_crystal_11 + g23 * O_crystal_12 + g33 * O_crystal_13;

% g_rot_symm_21 = g11 * O_crystal_21 + g21 * O_crystal_22 + g31 * O_crystal_23;
% g_rot_symm_22 = g12 * O_crystal_21 + g22 * O_crystal_22 + g32 * O_crystal_23;
g_rot_symm_23 = g13 * O_crystal_21 + g23 * O_crystal_22 + g33 * O_crystal_23;

g_rot_symm_31 = g11 * O_crystal_31 + g21 * O_crystal_32 + g31 * O_crystal_33;
g_rot_symm_32 = g12 * O_crystal_31 + g22 * O_crystal_32 + g32 * O_crystal_33;
g_rot_symm_33 = g13 * O_crystal_31 + g23 * O_crystal_32 + g33 * O_crystal_33;

PHI_symm = acos(g_rot_symm_33);
phi1_symm = atan2(g_rot_symm_31,-g_rot_symm_32);
phi2_symm = atan2(g_rot_symm_13,g_rot_symm_23);

tmp_ix = phi1_symm < 0;
phi1_symm(tmp_ix) = phi1_symm(tmp_ix) + 2*pi;

tmp_ix = phi2_symm < 0;
phi2_symm(tmp_ix) = phi2_symm(tmp_ix) + 2*pi;

for i = 1:num_points
    ix = (i-1)*N + 1 : i*N;
    phi1_f(ix) = phi1_symm(i,:);
    PHI_f(ix)  = PHI_symm(i,:);
    phi2_f(ix) = phi2_symm(i,:);
end

