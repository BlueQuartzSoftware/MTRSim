% unit_triangle_IPF_coords.m

% Daniel M. Sparkman
% 07/05/2017
% Research

% unit_triangle_IPF_coords.m

% this code takes as input the euler angles and returns as output the IPF
% coordinates in the unit triangle

function [X_fundamental,Y_fundamental] = unit_triangle_IPF_coords(phi1,PHI,phi2)

% hcp     fcc
crystal_system = 'hcp';


%% setup

% PHI = mod(PHI,pi/180*90);
% phi2 = mod(phi2,pi/180*60);


% set normal vector of specimen plane of interest in specimen coords
h = [0;
     0;
     1];

num_points = numel(phi1); 

%% develop symmetry operations

full_pole_figure = true;

if full_pole_figure

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
    

end

%% compute IPF coordinates
ix = phi2 == 0;
phi2(ix) = phi2(ix)+1e-15;

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

% rotate specimen plane normal vector into crystal coord representation
% % h_rot_1 = g11*h(1) + g12*h(2) + g13*h(3);
% % h_rot_2 = g21*h(1) + g22*h(2) + g23*h(3);
% % h_rot_3 = g31*h(1) + g32*h(2) + g33*h(3);

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
g_rot_symm_11 = g11 * O_crystal_11 + g21 * O_crystal_12 + g31 * O_crystal_13;
g_rot_symm_12 = g12 * O_crystal_11 + g22 * O_crystal_12 + g32 * O_crystal_13;
g_rot_symm_13 = g13 * O_crystal_11 + g23 * O_crystal_12 + g33 * O_crystal_13;

g_rot_symm_21 = g11 * O_crystal_21 + g21 * O_crystal_22 + g31 * O_crystal_23;
g_rot_symm_22 = g12 * O_crystal_21 + g22 * O_crystal_22 + g32 * O_crystal_23;
g_rot_symm_23 = g13 * O_crystal_21 + g23 * O_crystal_22 + g33 * O_crystal_23;

g_rot_symm_31 = g11 * O_crystal_31 + g21 * O_crystal_32 + g31 * O_crystal_33;
g_rot_symm_32 = g12 * O_crystal_31 + g22 * O_crystal_32 + g32 * O_crystal_33;
g_rot_symm_33 = g13 * O_crystal_31 + g23 * O_crystal_32 + g33 * O_crystal_33;

% rotate specimen plane normal vector into crystal coord representation (each should be num_pointsxN)
h_rot_symm_1 = g_rot_symm_11*h(1) + g_rot_symm_12*h(2) + g_rot_symm_13*h(3);
h_rot_symm_2 = g_rot_symm_21*h(1) + g_rot_symm_22*h(2) + g_rot_symm_23*h(3);
h_rot_symm_3 = g_rot_symm_31*h(1) + g_rot_symm_32*h(2) + g_rot_symm_33*h(3);

clear g_rot_symm_11 g_rot_symm_12 g_rot_symm_13 g_rot_symm_21 g_rot_symm_22 g_rot_symm_23 g_rot_symm_31 g_rot_symm_32 g_rot_symm_33 

% just examining half sphere because of symmetry over each side of the plane
ix = h_rot_symm_3 > 0;
h_rot_symm_1(ix) = h_rot_symm_1(ix)*-1;
h_rot_symm_2(ix) = h_rot_symm_2(ix)*-1;
h_rot_symm_3(ix) = h_rot_symm_3(ix)*-1;
        
Xall = h_rot_symm_1 ./ (1-h_rot_symm_3);
Yall = h_rot_symm_2 ./ (1-h_rot_symm_3);

% get indices for poles that land in unit stereographic triangle
ix1 = Xall >= 0;
ix2 = Yall <= Xall*tan(pi/180*30); 
ix3 = Yall >= 0;
ix_all = ix1 & ix2 & ix3;

X_fundamental = zeros(num_points,1);
Y_fundamental = zeros(num_points,1);
[val,ix] = max(ix_all,[],2);
for j = 1:num_points
    X_fundamental(j,1) = Xall(j,ix(j,1));
    Y_fundamental(j,1) = Yall(j,ix(j,1));
end

if 0
    [val,ix] = min(max(ix_all,[],2))
    figure(1); 
    plot(Xall(ix,:),Yall(ix,:),'b.'); 
    hold on; plot([0 1 cos(pi/180*30) 0],[0 0 sin(pi/180*30) 0],'k-'); hold off;
    axis([-1 1 -1 1]); axis square
end

