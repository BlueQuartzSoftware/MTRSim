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

% Eps-snap: match EbsdLib eu2om's rounding cleanup. Any matrix entry
% smaller than 1e-7 in magnitude is zeroed -- prevents tiny FP noise
% from polluting downstream matrix products and Euler extractions.
eps_snap = 1.0e-7;
g11(abs(g11) < eps_snap) = 0; g12(abs(g12) < eps_snap) = 0; g13(abs(g13) < eps_snap) = 0;
g21(abs(g21) < eps_snap) = 0; g22(abs(g22) < eps_snap) = 0; g23(abs(g23) < eps_snap) = 0;
g31(abs(g31) < eps_snap) = 0; g32(abs(g32) < eps_snap) = 0; g33(abs(g33) < eps_snap) = 0;

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

% Eps-snap: same cleanup applied to the symmetry-operator matrices so
% sin(0)/cos(pi/2) FP residue can't sneak through the product either.
O_crystal_11(abs(O_crystal_11) < eps_snap) = 0; O_crystal_12(abs(O_crystal_12) < eps_snap) = 0; O_crystal_13(abs(O_crystal_13) < eps_snap) = 0;
O_crystal_21(abs(O_crystal_21) < eps_snap) = 0; O_crystal_22(abs(O_crystal_22) < eps_snap) = 0; O_crystal_23(abs(O_crystal_23) < eps_snap) = 0;
O_crystal_31(abs(O_crystal_31) < eps_snap) = 0; O_crystal_32(abs(O_crystal_32) < eps_snap) = 0; O_crystal_33(abs(O_crystal_33) < eps_snap) = 0;


% build symmetry rotation matrices (each should be num_pointsxN)
g_rot_symm_11 = g11 * O_crystal_11 + g21 * O_crystal_12 + g31 * O_crystal_13;
g_rot_symm_12 = g12 * O_crystal_11 + g22 * O_crystal_12 + g32 * O_crystal_13;
g_rot_symm_13 = g13 * O_crystal_11 + g23 * O_crystal_12 + g33 * O_crystal_13;

% g_rot_symm_21 = g11 * O_crystal_21 + g21 * O_crystal_22 + g31 * O_crystal_23;
% g_rot_symm_22 = g12 * O_crystal_21 + g22 * O_crystal_22 + g32 * O_crystal_23;
g_rot_symm_23 = g13 * O_crystal_21 + g23 * O_crystal_22 + g33 * O_crystal_23;

g_rot_symm_31 = g11 * O_crystal_31 + g21 * O_crystal_32 + g31 * O_crystal_33;
g_rot_symm_32 = g12 * O_crystal_31 + g22 * O_crystal_32 + g32 * O_crystal_33;
g_rot_symm_33 = g13 * O_crystal_31 + g23 * O_crystal_32 + g33 * O_crystal_33;

% MATLAB om2eu -- matches EbsdLib::OrientationMatrix::toEuler convention.
% Three branches: non-degenerate (PHI strictly inside (0, pi)), PHI~=0
% (gimbal lock at the lower bound, recover c-axis rotation), PHI~=pi
% (gimbal lock at the upper bound, recover c-axis rotation).
%
% Without the degenerate branches the code silently throws away the
% c-axis rotation amount on gimbal-lock cases via atan2(0, 0) = 0.
% That mismatches both EbsdLib (which we trust as the MTEX-validated
% reference) and the C++ ComputeODFFilter implementation.

eps_close = 1.0e-6;
g33_abs = abs(g_rot_symm_33);
near_one = abs(g33_abs - 1.0) <= eps_close;
non_degen = ~near_one;
near_pos1 = near_one & (g_rot_symm_33 >= 0.0);
near_neg1 = near_one & (g_rot_symm_33 <  0.0);

PHI_symm  = zeros(size(g_rot_symm_33));
phi1_symm = zeros(size(g_rot_symm_33));
phi2_symm = zeros(size(g_rot_symm_33));

% --- Non-degenerate branch
PHI_symm(non_degen)  = acos(g_rot_symm_33(non_degen));
phi1_symm(non_degen) = atan2(g_rot_symm_31(non_degen), -g_rot_symm_32(non_degen));
phi2_symm(non_degen) = atan2(g_rot_symm_13(non_degen),  g_rot_symm_23(non_degen));

% --- PHI ~= 0 (gimbal lock at lower bound)
% C-axis rotation amount goes into phi1; PHI = 0; phi2 = 0.
phi1_symm(near_pos1) = atan2(g_rot_symm_12(near_pos1), g_rot_symm_11(near_pos1));
% PHI_symm and phi2_symm stay 0 from initialisation.

% --- PHI ~= pi (gimbal lock at upper bound)
% C-axis rotation amount goes into phi1 with the EbsdLib sign; PHI = pi; phi2 = 0.
phi1_symm(near_neg1) = -atan2(-g_rot_symm_12(near_neg1), g_rot_symm_11(near_neg1));
PHI_symm(near_neg1)  = pi;
% phi2_symm stays 0.

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

