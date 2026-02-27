function [PF] = convert_ODF_to_PF(ODFval, ODFbins, phi1_bins, PHI_bins, phi2_bins)

% convert_ODF_to_PF.m

% Daniel M. Sparkman
% 08/18/2013
% Research

% convert_ODF_to_PF.m

% this code takes an ODF for Euler angles (Bunge notation) and returns
% the intensity of the pole figure for each orientation bin

%% setup

% discretization of euler space
degree_spacing = 5;
num_bins = 360 / degree_spacing;

%% calc

% size of bins
dphi1 = phi1_bins(2) - phi1_bins(1);
dPHI = PHI_bins(2) - PHI_bins(1);
dphi2 = phi2_bins(2) - phi2_bins(1);

fprintf('calculating pole figure coordinates ...\n');
% fprintf('%%complete : 0%% ');

use_prebuilt_PF_coords = true;

if use_prebuilt_PF_coords
    load PF_coords_for_ODF PFval
    
    for m = 1:size(ODFval,1)
        % examining the m^th bin
        if ODFval(m,1) == 0
            PFval(m,3) = 0;
        else
            PFval(m,3) = ODFval(m,1) * 4*pi^2 / ( dphi1 * dphi2 * abs(cos(ODFbins(m,2)- dPHI/2) - cos(ODFbins(m,2)+ dPHI/2))) ;
        end
    end
else

PFval = zeros(size(ODFval,1),3);
for m = 1:size(ODFval,1)
    % examining the m^th bin
    % euler space position the bin is centered about
    EA_avg = ODFbins(m,1:3);
    
    % coordinates of pole in pole figure
    [X,Y] = pole_figure_coords(EA_avg(1), EA_avg(2), EA_avg(3));
    PFval(m,1) = X;
    PFval(m,2) = Y;
    
    if ODFval(m,1) == 0
        PFval(m,3) = 0;
    else
        PFval(m,3) = ODFval(m,1) * 4*pi^2 / ( dphi1 * dphi2 * abs(cos(ODFbins(m,2)- dPHI/2) - cos(ODFbins(m,2)+ dPHI/2))) ;
    end
end
end


%% find all of the euler space bins that fall within each pole figure bin
% (due to crystal symmetry)

use_prebuilt_index = true;

if use_prebuilt_index
%     load IX_ODF_to_PF IX_ODF_to_PF PF_indexed
    load PF_indexed IX_ODF_to_PF PF_indexed
    
    num_PF_bins = max(IX_ODF_to_PF);
    PFval_tmp = zeros(num_PF_bins,3);
    for j = 1:num_PF_bins
        bin_IX = IX_ODF_to_PF == j;
        PFval_tmp(j,:) = [PF_indexed(j,1:2) sum(PFval( bin_IX , 3))];
    end
    PFval = PFval_tmp;
else
PFval_tmp = [];
i = 0;
j = 0;
tol = 0.0001;
while i < size(PFval,1)
    % find indices of eulerspace bins that correspond to pole figure bins
    i = i+1;
    IXx = (PFval(:,1) >= PFval(i,1) -tol) & (PFval(:,1) <= PFval(i,1) +tol);
    IXy = (PFval(:,2) >= PFval(i,2) -tol) & (PFval(:,2) <= PFval(i,2) +tol);
    IX = IXx & IXy; % logical(IXx .* IXy);
    
    % increment pole figure intensity value
    j = j+1;
    PFval_tmp(j,:) = [PFval(i,1:2) sum(PFval( IX , 3))];
%     PFval_tmp(j,:) = [PFval(i,1:2) mean(PFval( IX , 3))];
    PFval(IX,:) = [];
    i = i-1;
end
PFval = PFval_tmp;
end

%% plot it
PF = struct('PFval',PFval,'num_bins',num_bins);
