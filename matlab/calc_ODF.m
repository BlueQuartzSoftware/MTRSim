function [ODFval, ODFbins, phi1_bins, PHI_bins, phi2_bins] = calc_ODF(phi1, PHI, phi2)

% calc_ODF.m

% Daniel M. Sparkman
% 08/13/2013
% Research

% calc_ODF.m

% this code takes a collection of Euler angles (Bunge notation) and returns
% the ODF value for each orientation bin

%% setup

% discretization of euler space
degree_spacing = 5;
radian_spacing = degree_spacing * pi / 180;
num_bins = 360 / degree_spacing;

% discretize euler space
phi1_bins = 0:2*pi/num_bins:2*pi;
PHI_bins  = 0:pi*2/num_bins:pi;
phi2_bins = 0:2*pi/num_bins:2*pi;

num_phi1_bins = length(phi1_bins);
num_PHI_bins = length(PHI_bins);
num_phi2_bins = length(phi2_bins);

num_total_bins = (num_phi1_bins-1)*(num_PHI_bins-1)*(num_phi2_bins-1);

% obtain all symmetric euler angles
if 1
    fprintf('calculating symmetric orientations...\n');
    [phi1_vec, PHI_vec, phi2_vec] = symmetric_euler_angles(phi1,PHI,phi2);
    n = size(phi1_vec,1);
else
    phi1_vec = phi1;
    PHI_vec = PHI;
    phi2_vec = phi2;
end

phi1_ix = fix(phi1_vec/radian_spacing)+1;
phi1_ix(phi1_ix == num_bins+1) = phi1_ix(phi1_ix == num_bins+1)-1;
PHI_ix  = fix(PHI_vec/radian_spacing)+1;
PHI_ix(PHI_ix == num_bins/2+1) = PHI_ix(PHI_ix == num_bins/2+1)-1;
phi2_ix = fix(phi2_vec/radian_spacing)+1;
phi2_ix(phi2_ix == num_bins+1) = phi2_ix(phi2_ix == num_bins+1)-1;

% load blank ODF
use_blank_ODF = true;

if use_blank_ODF
    load blank_ODF ODFbins ODFval phi1_bins PHI_bins phi2_bins
else
    m = 0;
    ODFbins = zeros((num_phi1_bins-1)*(num_PHI_bins-1)*(num_phi2_bins-1), 3);
    ODFval = zeros((num_phi1_bins-1)*(num_PHI_bins-1)*(num_phi2_bins-1), 1);
    for j = 2:num_phi1_bins
        for k = 2:num_PHI_bins
            for l = 2:num_phi2_bins
                % examining the m^th bin
                m = m+1;
                
                % euler space position the bin is centered about
                EA_avg = [ (phi1_bins(j-1)+phi1_bins(j))/2    (PHI_bins(k-1)+PHI_bins(k))/2    (phi2_bins(l-1)+phi2_bins(l))/2 ];
                ODFbins(m,1:3) = EA_avg;
            end
        end
    end
    save blank_ODF ODFbins ODFval phi1_bins PHI_bins phi2_bins
end

ODFval(:,1) = zeros(size(ODFval,1),1);
fprintf('indexing symmetric orientations...\n');
% indexing
ix_fg_symm = (phi1_ix-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (PHI_ix-1)*(num_phi2_bins-1) + phi2_ix;

% number of orientations
N = size(phi1_vec,1);

% set smoothing
smooth_ODF = true;
    
unique_fg_ix = unique(ix_fg_symm);
num_uniques = numel(unique_fg_ix);

%% calc

if N < num_total_bins || 1  % update orienation by orientation
    
    j = phi1_ix + 1;
    k =  PHI_ix + 1;
    l = phi2_ix + 1;
    
    jf = j - 1;
    kf = k - 1;
    lf = l - 1;
    
    jf_minus = jf -1;
    jf_minus(j == 2) = num_phi1_bins - 1;

    jf_plus = jf +1;
	jf_plus(j == num_phi1_bins) = 1;
        
    kf_minus = kf -1;
    kf_minus(k == 2) = num_PHI_bins - 1;
    kf_plus = kf +1;
    kf_plus(k == num_PHI_bins) = 1;
    
    lf_minus = lf -1;
    lf_minus(l == 2) = num_phi2_bins - 1;
    lf_plus = lf +1;
    lf_plus(l == num_phi2_bins) = 1;
    
    % update ODF for this bin
    center_P_factor = 0.332;
    
    % update ODF for bins neighboring this bin's faces
    face_P_factor = 0.448/6;
    ix_j_minus = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1) + (lf);
    ix_j_plus  = (jf_plus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1) + (lf);
    ix_k_minus = (jf-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf);
    ix_k_plus  = (jf-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1) + (lf);
    ix_l_minus = (jf-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1) + (lf_minus);
    ix_l_plus  = (jf-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1) + (lf_plus);
    ix_face_all = [ix_j_minus ix_j_plus ix_k_minus ix_k_plus ix_l_minus ix_l_plus];
    clear ix_j_minus ix_j_plus ix_k_minus ix_k_plus ix_l_minus ix_l_plus
    
    % update edge bins
    edge_P_factor = 0.16 / 12;
    ix_edge_1  = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1)       + (lf_minus);
    ix_edge_2  = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf);
    ix_edge_3  = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1)  + (lf);
    ix_edge_4  = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1)       + (lf_plus);
    ix_edge_5  = (jf-1)      *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf_minus);
    ix_edge_6  = (jf-1)      *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf_plus);
    ix_edge_7  = (jf-1)      *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1)  + (lf_plus);
    ix_edge_8  = (jf-1)      *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1)  + (lf_minus);
    ix_edge_9  = (jf_plus-1) *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1)       + (lf_minus);
    ix_edge_10 = (jf_plus-1) *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf);
    ix_edge_11 = (jf_plus-1) *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1)  + (lf);
    ix_edge_12 = (jf_plus-1) *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1)       + (lf_plus);
    ix_edge_all = [ix_edge_1 ix_edge_2 ix_edge_3 ix_edge_4 ix_edge_5 ix_edge_6 ix_edge_7 ix_edge_8 ix_edge_9 ix_edge_10 ix_edge_11 ix_edge_12];
    clear ix_edge_1 ix_edge_2 ix_edge_3 ix_edge_4 ix_edge_5 ix_edge_6 ix_edge_7 ix_edge_8 ix_edge_9 ix_edge_10 ix_edge_11 ix_edge_12
    
    % update corner bins
    corner_P_factor = 0.06 / 8;
    ix_corner_1 = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf_minus);
    ix_corner_2 = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf_plus);
    ix_corner_3 = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1) + (lf_plus);
    ix_corner_4 = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1) + (lf_minus);
    ix_corner_5 = (jf_plus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf_minus);
    ix_corner_6 = (jf_plus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf_plus);
    ix_corner_7 = (jf_plus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1) + (lf_plus);
    ix_corner_8 = (jf_plus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1) + (lf_minus);
    ix_corner_all = [ix_corner_1 ix_corner_2 ix_corner_3 ix_corner_4 ix_corner_5 ix_corner_6 ix_corner_7 ix_corner_8];
    clear ix_corner_1 ix_corner_2 ix_corner_3 ix_corner_4 ix_corner_5 ix_corner_6 ix_corner_7 ix_corner_8
    
    if smooth_ODF
        fprintf('calculating ODF ...\n');
        if 1
            itm = (1:num_total_bins)';
            count = hist(ix_fg_symm,itm)';
%             [val,ixc] = unique(ix_fg_symm);
            
            ODFval(itm,1) = ODFval(itm,1) + count / N * center_P_factor ;
            
            count = hist(ix_face_all(:),itm)';
            ODFval(itm,1) = ODFval(itm,1) + count / N * face_P_factor ;
            
            count = hist(ix_edge_all(:),itm)';
            ODFval(itm,1) = ODFval(itm,1) + count / N * edge_P_factor ;
            
            count = hist(ix_corner_all(:),itm)';
            ODFval(itm,1) = ODFval(itm,1) + count / N * corner_P_factor ;
            
        else
            fprintf('%% complete : 0%% ');
            disp_count = 0;
            disp_freq = ceil(N/19.99);
            for i = 1:N
                disp_count = disp_count+1;
                if disp_count == disp_freq
                    disp_count = 0;
                    fprintf('%g%% ',round(100*i/N));
                end
                
                ODFval(ix_fg_symm(i),1) = ODFval(ix_fg_symm(i),1) + 1/N * center_P_factor;
                for j = 1:size(ix_face_all,2)
                    ODFval(ix_face_all(i,j),1) = ODFval(ix_face_all(i,j),1) + 1/N * face_P_factor;
                end
                for j = 1:size(ix_edge_all,2)
                    ODFval(ix_edge_all(i,j),1)  = ODFval(ix_edge_all(i,j),1)  + 1/N *  edge_P_factor;
                end
                for j = 1:size(ix_corner_all,2)
                    ODFval(ix_corner_all(i,j),1) = ODFval(ix_corner_all(i,j),1) + 1/N *  corner_P_factor;
                end
            end
            fprintf('100%%\n');
        end
    else
        for i = 1:N
            ODFval(ix_fg_symm(i),1) = ODFval(ix_fg_symm(i),1) + 1/N;
        end
    end
    
    fprintf('\n');

    
elseif N < num_total_bins   % update orienation by orientation
    fg_ix = (phi1_ix-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (PHI_ix-1)*(num_phi2_bins-1) + phi2_ix;
    unique_fg_ix = unique(fg_ix);
    num_uniques = numel(unique_fg_ix);
    
    if use_blank_ODF
        load blank_ODF ODFbins ODFval phi1_bins PHI_bins phi2_bins
    else
        m = 0;
        ODFbins = zeros((num_phi1_bins-1)*(num_PHI_bins-1)*(num_phi2_bins-1), 3);
        ODFval = zeros((num_phi1_bins-1)*(num_PHI_bins-1)*(num_phi2_bins-1), 1);
        for j = 2:num_phi1_bins
            for k = 2:num_PHI_bins
                for l = 2:num_phi2_bins
                    % examining the m^th bin
                    m = m+1;
                    
                    % euler space position the bin is centered about
                    EA_avg = [ (phi1_bins(j-1)+phi1_bins(j))/2    (PHI_bins(k-1)+PHI_bins(k))/2    (phi2_bins(l-1)+phi2_bins(l))/2 ];
                    ODFbins(m,1:3) = EA_avg;
                end
            end
        end
        save blank_ODF ODFbins ODFval phi1_bins PHI_bins phi2_bins
    end
    
    fprintf('calculating ODF ...\n');
    fprintf('%% complete : 0%% ');
    
    disp_count = 0;
    disp_freq = ceil(num_uniques/9.99);
    
    for m = 1 : num_uniques
        disp_count = disp_count+1;
        if disp_count == disp_freq
            disp_count = 0;
            fprintf('%g%% ',round(100*m/num_uniques));
        end
        
        ix = unique_fg_ix(m); % index of bin of interest
        IX = fg_ix == ix; % indices of orientations that fall in this bin
        sumIX = sum(IX); % total number of orientations that fall in this bin
        
        [val,IX_tmp] = max(IX);
        phi1_ix_single = phi1_ix(IX_tmp);
        PHI_ix_single  = PHI_ix(IX_tmp);
        phi2_ix_single = phi2_ix(IX_tmp);
        
        % smooth the ODF by incrementing neighboring bins as well
        if smooth_ODF
            if 0
                tmp = fix( ix/((num_PHI_bins-1)*(num_phi2_bins-1)) );
                if tmp == num_phi1_bins - 1
                    j = tmp+2;
                    remainder = ((num_PHI_bins-1)*(num_phi2_bins-1));
                else
                    j = fix( ix/((num_PHI_bins-1)*(num_phi2_bins-1)) ) + 2;
                    remainder = rem( ix, ((num_PHI_bins-1)*(num_phi2_bins-1)) );
                end
                k = fix( remainder/(num_phi2_bins-1) ) + 2;
                l = rem(remainder, (num_phi2_bins-1)) + 1;
            else
                j = phi1_ix_single + 1;
                k = PHI_ix_single  + 1;
                l = phi2_ix_single + 1;
            end

            jf = j - 1;
            kf = k - 1;
            lf = l - 1;
            m_check = ((j-1)-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (k-2)*(num_phi2_bins-1) + (l-1);
            
            if j == 2
                jf_minus = num_phi1_bins - 1;
            else
                jf_minus = jf -1;
            end
            if j == num_phi1_bins
                jf_plus = 1;
            else
                jf_plus = jf +1;
            end
            
            if k == 2
                kf_minus = num_PHI_bins - 1;
            else
                kf_minus = kf -1;
            end
            if k == num_PHI_bins
                kf_plus = 1;
            else
                kf_plus = kf +1;
            end
            
            if l == 2
                lf_minus = num_phi2_bins - 1;
            else
                lf_minus = lf -1;
            end
            if l == num_phi2_bins
                lf_plus = 1;
            else
                lf_plus = lf +1;
            end
            
            % update ODF for this bin
            center_P_factor = 0.332;
            ODFval(ix,1) = ODFval(ix,1) + sumIX/N * center_P_factor;
            
            % update ODF for bins neighboring this bin's faces
            face_P_factor = 0.448/6;
            ix_j_minus = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1) + (lf);
            ix_j_plus  = (jf_plus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1) + (lf);
            ix_k_minus = (jf-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf);
            ix_k_plus  = (jf-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1) + (lf);
            ix_l_minus = (jf-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1) + (lf_minus);
            ix_l_plus  = (jf-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1) + (lf_plus);
            ODFval(ix_j_minus,1) = ODFval(ix_j_minus,1) + sumIX/N * face_P_factor;
            ODFval(ix_j_plus,1)  = ODFval(ix_j_plus,1)  + sumIX/N * face_P_factor;
            ODFval(ix_k_minus,1) = ODFval(ix_k_minus,1) + sumIX/N * face_P_factor;
            ODFval(ix_k_plus,1)  = ODFval(ix_k_plus,1)  + sumIX/N * face_P_factor;
            ODFval(ix_l_minus,1) = ODFval(ix_l_minus,1) + sumIX/N * face_P_factor;
            ODFval(ix_l_plus,1)  = ODFval(ix_l_plus,1)  + sumIX/N * face_P_factor;
            
            % update edge bins
            edge_P_factor = 0.16 / 12;
            ix_edge_1  = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1)       + (lf_minus);
            ix_edge_2  = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf);
            ix_edge_3  = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1)  + (lf);
            ix_edge_4  = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1)       + (lf_plus);
            ix_edge_5  = (jf-1)      *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf_minus);
            ix_edge_6  = (jf-1)      *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf_plus);
            ix_edge_7  = (jf-1)      *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1)  + (lf_plus);
            ix_edge_8  = (jf-1)      *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1)  + (lf_minus);
            ix_edge_9  = (jf_plus-1) *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1)       + (lf_minus);
            ix_edge_10 = (jf_plus-1) *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf);
            ix_edge_11 = (jf_plus-1) *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1)  + (lf);
            ix_edge_12 = (jf_plus-1) *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1)       + (lf_plus);
            ODFval(ix_edge_1,1)  = ODFval(ix_edge_1,1)  + sumIX/N *  edge_P_factor;
            ODFval(ix_edge_2,1)  = ODFval(ix_edge_2,1)  + sumIX/N *  edge_P_factor;
            ODFval(ix_edge_3,1)  = ODFval(ix_edge_3,1)  + sumIX/N *  edge_P_factor;
            ODFval(ix_edge_4,1)  = ODFval(ix_edge_4,1)  + sumIX/N *  edge_P_factor;
            ODFval(ix_edge_5,1)  = ODFval(ix_edge_5,1)  + sumIX/N *  edge_P_factor;
            ODFval(ix_edge_6,1)  = ODFval(ix_edge_6,1)  + sumIX/N *  edge_P_factor;
            ODFval(ix_edge_7,1)  = ODFval(ix_edge_7,1)  + sumIX/N *  edge_P_factor;
            ODFval(ix_edge_8,1)  = ODFval(ix_edge_8,1)  + sumIX/N *  edge_P_factor;
            ODFval(ix_edge_9,1)  = ODFval(ix_edge_9,1)  + sumIX/N *  edge_P_factor;
            ODFval(ix_edge_10,1) = ODFval(ix_edge_10,1) + sumIX/N *  edge_P_factor;
            ODFval(ix_edge_11,1) = ODFval(ix_edge_11,1) + sumIX/N *  edge_P_factor;
            ODFval(ix_edge_12,1) = ODFval(ix_edge_12,1) + sumIX/N *  edge_P_factor;
            
            % update corner bins
            corner_P_factor = 0.06 / 8;
            ix_corner_1 = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf_minus);
            ix_corner_2 = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf_plus);
            ix_corner_3 = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1) + (lf_plus);
            ix_corner_4 = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1) + (lf_minus);
            ix_corner_5 = (jf_plus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf_minus);
            ix_corner_6 = (jf_plus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf_plus);
            ix_corner_7 = (jf_plus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1) + (lf_plus);
            ix_corner_8 = (jf_plus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1) + (lf_minus);
            ODFval(ix_corner_1,1) = ODFval(ix_corner_1,1) + sumIX/N *  corner_P_factor;
            ODFval(ix_corner_2,1) = ODFval(ix_corner_2,1) + sumIX/N *  corner_P_factor;
            ODFval(ix_corner_3,1) = ODFval(ix_corner_3,1) + sumIX/N *  corner_P_factor;
            ODFval(ix_corner_4,1) = ODFval(ix_corner_4,1) + sumIX/N *  corner_P_factor;
            ODFval(ix_corner_5,1) = ODFval(ix_corner_5,1) + sumIX/N *  corner_P_factor;
            ODFval(ix_corner_6,1) = ODFval(ix_corner_6,1) + sumIX/N *  corner_P_factor;
            ODFval(ix_corner_7,1) = ODFval(ix_corner_7,1) + sumIX/N *  corner_P_factor;
            ODFval(ix_corner_8,1) = ODFval(ix_corner_8,1) + sumIX/N *  corner_P_factor;
            
        else
            ODFval(ix,1) = sumIX/N ;
        end
        
    end
    fprintf('100%%\n');
else  % update bin by bin
    
% size of bins
% dphi1 = phi1_bins(2) - phi1_bins(1);
% dPHI = PHI_bins(2) - PHI_bins(1);
% dphi2 = phi2_bins(2) - phi2_bins(1);
fprintf('calculating ODF ...\n');
fprintf('%% complete : 0%% ');
m = 0;
ODFval = zeros((num_phi1_bins-1)*(num_PHI_bins-1)*(num_phi2_bins-1), 1);
ODFbins = zeros((num_phi1_bins-1)*(num_PHI_bins-1)*(num_phi2_bins-1), 3);
disp_count = 0;
disp_freq = ceil(length(phi1_bins)/19.99);
for j = 2:num_phi1_bins
    disp_count = disp_count+1;
    if disp_count == disp_freq
        disp_count = 0;
        fprintf('%g%% ',round(100*j/length(phi1_bins)));
    end
    
    for k = 2:num_PHI_bins
        for l = 2:num_phi2_bins
            % examining the m^th bin
            m = m+1;
            
            % euler space position the bin is centered about
            EA_avg = [ (phi1_bins(j-1)+phi1_bins(j))/2    (PHI_bins(k-1)+PHI_bins(k))/2    (phi2_bins(l-1)+phi2_bins(l))/2 ];
            
            ODFbins(m,1:3) = EA_avg;
            
            % indices of orientations that fall within the m^th bin
%             IX1 = (phi1_vec > phi1_bins(j-1)) & (phi1_vec < phi1_bins(j));
%             IX2 = (PHI_vec  > PHI_bins(k-1))  & (PHI_vec  < PHI_bins(k));
%             IX3 = (phi2_vec > phi2_bins(l-1)) & (phi2_vec < phi2_bins(l));

%             sumIX = 0;
%             IX1 = phi1_ix == j-1;
%             if any(IX1)
%                 IX2 = PHI_ix == k-1;
%                 if any(IX2)
%                     IX3 = phi2_ix == l-1;
%                     if any(IX3)
%                         IX  = IX1 .* IX2 .* IX3;
%                         % number of orientations in the m^th bin
%                         sumIX = sum(IX);
%                     end
%                 end
%             end
            
            IX = (phi1_ix == j-1) & (PHI_ix == k-1) & (phi2_ix == l-1);
            sumIX = sum(IX);
            
            % only update the ODF for this bin if it contains observations 
            if sumIX > 0
                % comment out the normalization to speed up calculation
%                 PFval(m,3) = sumIX/N * length(phi1_bins) * length(PHI_bins) * length(phi2_bins) / (sin((PHI_bins(k) + PHI_bins(k-1))/2) * dphi1 * dPHI * dphi2) ;
%                 PFval(m,3) = sumIX/N / ((1/(8*pi^2)) * sin((PHI_bins(k) + PHI_bins(k-1))/2) * dphi1 * dPHI * dphi2 ) ;
%                 ODFval(m,1) = sumIX/N * 4*pi^2 / ( dphi1 * dphi2 * (cos(PHI_bins(k-1)) - cos(PHI_bins(k)))) ;
                
                % smooth the ODF by incrementing neighboring bins as well
                if smooth_ODF
                    jf = j - 1;
                    kf = k - 1;
                    lf = l - 1;
                    m_check = ((j-1)-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (k-2)*(num_phi2_bins-1) + (l-1);
                    
                    if j == 2
                        jf_minus = num_phi1_bins - 1;
                    else
                        jf_minus = jf -1;
                    end
                    if j == num_phi1_bins
                        jf_plus = 1;
                    else
                        jf_plus = jf +1;
                    end
                    
                    if k == 2
                        kf_minus = num_PHI_bins - 1;
                    else
                        kf_minus = kf -1;
                    end
                    if k == num_PHI_bins
                        kf_plus = 1;
                    else
                        kf_plus = kf +1;
                    end
                    
                    if l == 2
                        lf_minus = num_phi2_bins - 1;
                    else
                        lf_minus = lf -1;
                    end
                    if l == num_phi2_bins
                        lf_plus = 1;
                    else
                        lf_plus = lf +1;
                    end
                    
                    % update ODF for this bin
                    center_P_factor = 0.332;
                    ODFval(m,1) = ODFval(m,1) + sumIX/N * center_P_factor;
                    
                    % update ODF for bins neighboring this bin's faces
                    face_P_factor = 0.448/6;
                    ix_j_minus = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1) + (lf);
                    ix_j_plus  = (jf_plus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1) + (lf);
                    ix_k_minus = (jf-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf);
                    ix_k_plus  = (jf-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1) + (lf);
                    ix_l_minus = (jf-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1) + (lf_minus);
                    ix_l_plus  = (jf-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1) + (lf_plus);
                    ODFval(ix_j_minus,1) = ODFval(ix_j_minus,1) + sumIX/N * face_P_factor;
                    ODFval(ix_j_plus,1)  = ODFval(ix_j_plus,1)  + sumIX/N * face_P_factor;
                    ODFval(ix_k_minus,1) = ODFval(ix_k_minus,1) + sumIX/N * face_P_factor;
                    ODFval(ix_k_plus,1)  = ODFval(ix_k_plus,1)  + sumIX/N * face_P_factor;
                    ODFval(ix_l_minus,1) = ODFval(ix_l_minus,1) + sumIX/N * face_P_factor;
                    ODFval(ix_l_plus,1)  = ODFval(ix_l_plus,1)  + sumIX/N * face_P_factor;
                    
                    % update edge bins
                    edge_P_factor = 0.16 / 12;
                    ix_edge_1  = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1)       + (lf_minus);
                    ix_edge_2  = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf);
                    ix_edge_3  = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1)  + (lf);
                    ix_edge_4  = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1)       + (lf_plus);
                    ix_edge_5  = (jf-1)      *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf_minus);
                    ix_edge_6  = (jf-1)      *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf_plus);
                    ix_edge_7  = (jf-1)      *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1)  + (lf_plus);
                    ix_edge_8  = (jf-1)      *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1)  + (lf_minus);
                    ix_edge_9  = (jf_plus-1) *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1)       + (lf_minus);
                    ix_edge_10 = (jf_plus-1) *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf);
                    ix_edge_11 = (jf_plus-1) *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1)  + (lf);
                    ix_edge_12 = (jf_plus-1) *(num_PHI_bins-1)*(num_phi2_bins-1) + (kf-1)*(num_phi2_bins-1)       + (lf_plus);
                    ODFval(ix_edge_1,1)  = ODFval(ix_edge_1,1)  + sumIX/N *  edge_P_factor;
                    ODFval(ix_edge_2,1)  = ODFval(ix_edge_2,1)  + sumIX/N *  edge_P_factor;
                    ODFval(ix_edge_3,1)  = ODFval(ix_edge_3,1)  + sumIX/N *  edge_P_factor;
                    ODFval(ix_edge_4,1)  = ODFval(ix_edge_4,1)  + sumIX/N *  edge_P_factor;
                    ODFval(ix_edge_5,1)  = ODFval(ix_edge_5,1)  + sumIX/N *  edge_P_factor;
                    ODFval(ix_edge_6,1)  = ODFval(ix_edge_6,1)  + sumIX/N *  edge_P_factor;
                    ODFval(ix_edge_7,1)  = ODFval(ix_edge_7,1)  + sumIX/N *  edge_P_factor;
                    ODFval(ix_edge_8,1)  = ODFval(ix_edge_8,1)  + sumIX/N *  edge_P_factor;
                    ODFval(ix_edge_9,1)  = ODFval(ix_edge_9,1)  + sumIX/N *  edge_P_factor;
                    ODFval(ix_edge_10,1) = ODFval(ix_edge_10,1) + sumIX/N *  edge_P_factor;
                    ODFval(ix_edge_11,1) = ODFval(ix_edge_11,1) + sumIX/N *  edge_P_factor;
                    ODFval(ix_edge_12,1) = ODFval(ix_edge_12,1) + sumIX/N *  edge_P_factor;
                    
                    % update corner bins
                    corner_P_factor = 0.06 / 8;
                    ix_corner_1 = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf_minus);
                    ix_corner_2 = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf_plus);
                    ix_corner_3 = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1) + (lf_plus);
                    ix_corner_4 = (jf_minus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1) + (lf_minus);
                    ix_corner_5 = (jf_plus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf_minus);
                    ix_corner_6 = (jf_plus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_minus-1)*(num_phi2_bins-1) + (lf_plus);
                    ix_corner_7 = (jf_plus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1) + (lf_plus);
                    ix_corner_8 = (jf_plus-1)*(num_PHI_bins-1)*(num_phi2_bins-1) + (kf_plus-1)*(num_phi2_bins-1) + (lf_minus);
                    ODFval(ix_corner_1,1) = ODFval(ix_corner_1,1) + sumIX/N *  corner_P_factor;
                    ODFval(ix_corner_2,1) = ODFval(ix_corner_2,1) + sumIX/N *  corner_P_factor;
                    ODFval(ix_corner_3,1) = ODFval(ix_corner_3,1) + sumIX/N *  corner_P_factor;
                    ODFval(ix_corner_4,1) = ODFval(ix_corner_4,1) + sumIX/N *  corner_P_factor;
                    ODFval(ix_corner_5,1) = ODFval(ix_corner_5,1) + sumIX/N *  corner_P_factor;
                    ODFval(ix_corner_6,1) = ODFval(ix_corner_6,1) + sumIX/N *  corner_P_factor;
                    ODFval(ix_corner_7,1) = ODFval(ix_corner_7,1) + sumIX/N *  corner_P_factor;
                    ODFval(ix_corner_8,1) = ODFval(ix_corner_8,1) + sumIX/N *  corner_P_factor;
                    
                else
                    ODFval(m,1) = sumIX/N ;
                end
            end
            
%             phi1_vec(logical(IX)) = [];
%             PHI_vec(logical(IX)) = [];
%             phi2_vec(logical(IX)) = [];
        end
    end
end
fprintf('\n');
% ODFval(:,1) = ODFval(:,1) / sum(ODFval(:,1));
end

%% plot it

if 0
    % plot it!
    x = ODFval(:,1);
    y = ODFval(:,2);
    
%     figure(1)
    [TRI] = delaunay(x,y);
%     voronoi(x,y)
    [v,c]=voronoin([x y]);
    
    cm = colormap('jet');
    cmin = 0;
    cmax = max(ODFval(:,1));
    caxis([cmin cmax])
    
    % plot the ODF values according to colorbar
    for i = 1:length(c)
        
        ci1 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,1),ODFval(i,1),'linear');
        ci2 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,2),ODFval(i,1),'linear');
        ci3 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,3),ODFval(i,1),'linear');
        
        tcolor(1,i,1:3) = [ci1 ci2 ci3];  % pole figure intensity
        
        if use_voronoi
            patch( v(c{i},1) , v(c{i},2) , tcolor(1,i,1:3), 'EdgeAlpha',0 );
        else
            hold on;
            plot(x(i),y(i),'o','MarkerFaceColor',tcolor(1,i,:),'MarkerEdgeColor',tcolor(1,i,:));
            hold off;
        end
    end
    
    colorbar
    caxis([cmin cmax])
    colorbar off
    
    set(gcf,'color','w');
    axis off;
end
