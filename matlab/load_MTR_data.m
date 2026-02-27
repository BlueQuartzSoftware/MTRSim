% load_MTR_data.m

% Daniel M. Sparkman
% 04/10/2018
% Research

% load_MTR_data.m

% this code loads the MTR data and makes an IPF map

clc
clear
close all

%% setup
curr_dir = pwd;

path_oi = 'C:\Users\'
cd(path_oi);

BoolMTR = csvread('BoolMTR.csv');   BoolMTR = logical(BoolMTR);
EulerAngles = csvread('EulerAngles.csv');
ParentIds = csvread('ParentIds.csv');
X_Position = csvread('X_Position.csv');
Y_Position = csvread('Y_Position.csv');

s_sim = [X_Position, Y_Position];
phi1_vec = EulerAngles(:,1);
PHI_vec = EulerAngles(:,2);
phi2_vec = EulerAngles(:,3);

%% make IPF map

a = view_IPF_map(s_sim, phi1_vec, PHI_vec, phi2_vec);
% imwrite(a,'sim_IPF_map.jpg','jpg');

%% Bool Maps

N = size(BoolMTR,1);
x_vec = s_sim(:,1);
y_vec = s_sim(:,2);

xv = unique(x_vec);
yv = unique(y_vec);

dx = xv(2)-xv(1);
dy = yv(2)-yv(1);

nx = 1e4;    nx = min([numel(xv), nx]);
ny = 1e4;    ny = min([numel(yv), ny]);

xmin = min(x_vec);
xmax = max(x_vec);
ymin = min(y_vec);
ymax = max(y_vec);

xlen = xmax-xmin;
ylen = ymax-ymin;

nx = round(xlen/dx) + 1;
ny = round(ylen/dy) + 1;
CI_mask = zeros(ny,nx);
for k = 1:N
    xt = x_vec(k,1);
    yt = y_vec(k,1);
    i = round(xt/dx);
    if xmin == 0
        i = i + 1;
    end
    j = round(yt/dy);
    if ymin == 0
        j = j + 1;
    end
    
    CI_mask(j, i) = BoolMTR(k);
end
% image([min(x_vec) max(x_vec)], [min(y_vec) max(y_vec)], CI_mask);
figure(88);
imshow(CI_mask);

%% Assignment Maps
if 1
N = size(BoolMTR,1);
x_vec = s_sim(:,1);
y_vec = s_sim(:,2);

xv = unique(x_vec);
yv = unique(y_vec);

dx = xv(2)-xv(1);
dy = yv(2)-yv(1);

nx = 1e4;    nx = min([numel(xv), nx]);
ny = 1e4;    ny = min([numel(yv), ny]);

xmin = min(x_vec);
xmax = max(x_vec);
ymin = min(y_vec);
ymax = max(y_vec);

xlen = xmax-xmin;
ylen = ymax-ymin;

nx = round(xlen/dx) + 1;
ny = round(ylen/dy) + 1;
ParentIDs_map = zeros(ny,nx);
for k = 1:N
    xt = x_vec(k,1);
    yt = y_vec(k,1);
    i = round(xt/dx);
    if xmin == 0
        i = i + 1;
    end
    j = round(yt/dy);
    if ymin == 0
        j = j + 1;
    end
    
    ParentIDs_map(j, i) = ParentIds(k);
end
% image([min(x_vec) max(x_vec)], [min(y_vec) max(y_vec)], CI_mask);
figure(89);
imshow(ParentIDs_map,[0,122899]);
end

%% find biggest MTR

numIDs = size(unique(ParentIds),1);

count_data = zeros(2,numIDs);
IDs = sort(unique(ParentIds))';
[th,IDs] = hist(ParentIds,IDs);
count_data(1,:) = IDs;
count_data(2,:) = th;
sorted_count_data = flipud(sortrows(count_data',2))';

%% plot MTRs of interest

ParentIDoi = sorted_count_data(1,1);

IXmap = ParentIDs_map == ParentIDoi;

numMTRs = 15;
MTRixs = [sorted_count_data(1,1), sorted_count_data(1,3:numMTRs+1)];

% MTRixs = MTRixs(15);
numMTRs = numel(MTRixs);
mappedIXmap = (numMTRs+1)*ones(size(IXmap));
for i = 1:numMTRs
    mappedIXmap(ParentIDs_map == MTRixs(i)) = i;
end
cmap = [1,0,0;
        0,1,0;
        0,0,1;
        0.8,0.8,0;
        0.3,1.0,0.3;
        rand(256-5,3)];
cmap(numMTRs+1,:) = [0,0,0];
cmap(numMTRs+2:end,:) = [];

a = zeros(ny,nx,3);
for i = 1:ny
    for j = 1:nx
        c = cmap( mappedIXmap(i,j), :);
        a(i,j,:) = c;
    end
end

figure(90);
% imshow(double(IXmap));
% imshow(uint8(mappedIXmap*255),cmap);
imshow(a);

% ODF
% IX = (ParentIds~=112720) & ~BoolMTR;
IX = ~BoolMTR;

IX = false(size(IX));
for i = 1:numMTRs
    ParentIDoi = MTRixs(i);
    IX = IX | (ParentIds == ParentIDoi);
end

[ODFval, ODFbins, phi1_bins, PHI_bins, phi2_bins] = calc_ODF(phi1_vec(IX), PHI_vec(IX), phi2_vec(IX));
ODF_macro = struct('ODFval',ODFval,'ODFbins',ODFbins,'phi1_bins',phi1_bins,'PHI_bins',PHI_bins,'phi2_bins',phi2_bins);
[PFval] = convert_ODF_to_PF(ODFval, ODFbins, phi1_bins, PHI_bins, phi2_bins);

figure(9);
plot_PF(PFval);
colorbar
axis square
title('macrotexture');

%% build ODF for all of the MTRs

IX = ~BoolMTR;

IX = false(size(IX));

ODF_mtx = zeros(size(ODFval,1),numMTRs);
for i = 1:numMTRs
    ParentIDoi = MTRixs(i);
    IX = (ParentIds == ParentIDoi);
    [ODFval, ODFbins, phi1_bins, PHI_bins, phi2_bins] = calc_ODF(phi1_vec(IX), PHI_vec(IX), phi2_vec(IX));
    ODF_mtx(:,i) = ODFval;
end

%% assign class to each MTR
MTR_class_ID = zeros(1,numMTRs);
num_MTR_classes = 3;
MTR_class_ODF = zeros(size(ODFval,1),num_MTR_classes);
MTR_class_ODF(:,1) = ODF_mtx(:,1);
MTR_class_ODF(:,2) = ODF_mtx(:,2);
MTR_class_ODF(:,3) = ODF_mtx(:,8);
% MTR_class_ODF(:,4) = ODF_mtx(:,8);
% MTR_class_ODF(:,5) = ODF_mtx(:,8);

MTR_class_assignment = zeros(ny,nx);
for i = 1:numMTRs
    val = inf;
    for j = 1:num_MTR_classes
        ODF_oi = ODF_mtx(:,i);
        ODF_class_test = MTR_class_ODF(:,j);
        theta = (ODF_oi' * ODF_class_test)/norm(ODF_oi)/norm(ODF_class_test);
        if theta < val
            val = theta;
            MTR_class_ID(1,i) = j;
        end
    end
    ParentIDoi = MTRixs(i);
    IX = (ParentIDs_map == ParentIDoi);
    MTR_class_assignment(IX) = MTR_class_ID(1,i);
end
figure(11);
imagesc(MTR_class_assignment/num_MTR_classes*255);

%%
i = 1;
[PFval] = convert_ODF_to_PF(ODF_mtx(:,i), ODFbins, phi1_bins, PHI_bins, phi2_bins);

figure(9);
plot_PF(PFval);
colorbar
axis square
title('macrotexture');
