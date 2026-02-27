% simulate_MTRs.m

% Daniel M. Sparkman
% 05/31/2017
% Research

% simulate_MTRs.m

% this code simulates MTRs using the plurigaussian model and FMM ODF.

clc
clear
close all

%% setup
tic;

save_new = true;
plot_on = true;

% load PGRF_parameters P theta_list nugvar s

% % set simulation volume parameters
%-------------------------
% length in each spatial direction  [mm]
xlen = 1.5 * 25.4; 
ylen = 0.5 * 25.4; 
zlen = 0.00 * 25.4;

% voxel spacing in each direction [mm]
dx = 0.02;
dy = 0.02;
dz = 0.02;

% % MTR spatial parameters
%-------------------------
% volume_fractions:   volume fraction of each MTR-/component-ODF in fullODF
volume_fractions = [0.30, 0.35, 0.35];

% theta_list:  correlation length parameters (columns: spatial direction)
%                                             rows:    latent GRF
%              note:  increasing one of these values will generally
%              increase the size of the MTRs in that spatial direction
theta_list = [0.10, 0.45, 0.10
              0.08, 0.37, 0.08]*1;

% nugget variance:   increases the noise in the simulation(may not be used)
nugvar = [0.67, 0.71, 0.72];

% %
%-------------------------
nx = round(xlen/dx);
ny = round(ylen/dy);
nz = max([round(zlen/dz),1]);
N = nx*ny*nz;
P = volume_fractions;

fprintf('simulation field : \n');
% fprintf(' x length: %4.3f cm       nx: %g       dx: %g um\n', xlen/1e4, nx, dx);
% fprintf(' y length: %4.3f cm       ny: %g       dx: %g um\n', ylen/1e4, ny, dy);
fprintf(' x length: %4.3f mm       nx: %g       dx: %g mm\n', xlen, nx, dx);
fprintf(' y length: %4.3f mm       ny: %g       dy: %g mm\n', ylen, ny, dy);
fprintf(' z length: %4.3f mm       nz: %g       dz: %g mm\n', zlen, nz, dz);
fprintf(' N : %g\n', N);
fprintf(' \n');

%% simulate plurigaussian random field

if save_new
    fprintf('simulating random draw of categorical random field...\n');
    s_sim = zeros(N,3);
    k = 0;
    for zix = 1:nz
        z_tmp = zix*dz;
        for j = 1:nx
            x_tmp = j*dx;
            for i = 1:ny
                y_tmp = i*dy;
                
                k = k+1;
                s_sim(k,:) = [ x_tmp y_tmp z_tmp];
            end
        end
    end
    
    [MTR_IX_sim, Y_sim] = PGRF_simulation(P, theta_list, nugvar, dx, dy, dz, nx, ny, nz);
    N = size(s_sim,1);
    volume_element_side_length_sim = max( max(s_sim,[],1) - min(s_sim,[],1) );
    fprintf(' complete.\n');

    fprintf(' writing to file...\n');
    save simulated_PGRF s_sim MTR_IX_sim Y_sim
    fprintf('     complete.\n');
else
    load simulated_PGRF s_sim MTR_IX_sim Y_sim
end

%%% plot underlying gaussians and assignment map
num_components = numel(P); 

% coloring for plot
custom_map = [0.1 0.1 0.1;
    1 0.9 0;
    1 0 1;
    0 1 0.1;
    0 0 1;
    1 1 0;
    0.25 0.25 0.75;
    0.75 0.25 0.75];
custom_map = [custom_map;
              rand(num_components-size(custom_map,1),3)];
% custom_map = [0.1 0.1 0.1;
%     1 0.9 0;
%     1 0 1];
if num_components < size(custom_map,1)
    custom_map(num_components+1:end,:) = [];
end

nx_sim = nx;
ny_sim = ny;
nz_sim = nz;
if 0  % plot underlying gaussians
    bounds = max([abs(max(Y_sim)) abs(min(Y_sim))]);  bounds = 3.5;
    for j = 1:size(Y_sim,2)
        figure(57+j)
        cm = colormap('jet');
        cmin = -bounds;
        cmax = bounds;
        caxis([cmin cmax])
        
        nx_sim = nx;
        ny_sim = ny;
        im = zeros(ny_sim,nx_sim,3);
        k = 0;
        for b = 1:nx_sim
            for a = 1:ny_sim
                k = k+1;
                ci1 = interp1(cmin:(cmax-cmin)/(size(cm,1)-1):cmax,cm(:,1)',Y_sim(k,j),'linear');
                ci2 = interp1(cmin:(cmax-cmin)/(size(cm,1)-1):cmax,cm(:,2)',Y_sim(k,j),'linear');
                ci3 = interp1(cmin:(cmax-cmin)/(size(cm,1)-1):cmax,cm(:,3)',Y_sim(k,j),'linear');
                tcolor = [ci1 ci2 ci3];  % microtexture assignment
                im(ny_sim+1-a,b,:) = tcolor;
            end
        end
        image(xlen*[0 1],ylen*[0 1],im);
        axis equal
        axis([xlen*[0 1], ylen*[0 1]])
        
        colorbar
        caxis([cmin cmax])
        colorbar off
        
        set(gcf,'color','w');
        axis off;
        hold on;
        plot(xlen*[0 1 1 0 0], ylen*[0 0 1 1 0],'k-');
        hold off;
        title(['Y_',num2str(j),'(s)'],'FontName','Times New Roman','FontSize',14,'FontWeight','b');
    end
end

%%% plot simulated assignment map

a = zeros(ny_sim,nx_sim,3);
k = 0;
zoi = 1;
for zix = 1:nz_sim
    for j = 1:nx_sim
        for i = 1:ny_sim
                k = k+1;
            if zix == zoi
                a(ny_sim+1-i,j,:) = custom_map(MTR_IX_sim(k),:);
            end
        end
    end
end
figure(105);
image([0 xlen],[0 ylen], a);
set(gcf,'color','w');
axis equal
cm = colormap(custom_map);
cmin = min(MTR_IX_sim(:));
cmax = max(MTR_IX_sim(:));
colorbar('Ticks',unique(MTR_IX_sim))
caxis([cmin cmax])

%% simulate orientations

if save_new

    load uniformODF
    load simulation_ODF ODF_best
    
    % normalize to integrate to 1
    for j = 1:numel(ODF_best)
        ODF_best(j).ODFval = ODF_best(j).ODFval/sum(ODF_best(j).ODFval); 
    end
    
    if 0
        ODF_best = ODF; N = 10000;
        tmp = ODF_best(5);
        for i = 1:numel(ODF)-1
            ODF_best(5+1-i) = ODF_best(5-i);
        end
        ODF_best(1) = tmp;
    end
    
    fprintf('\nsimulating orientations...(%g)\n',N);
    phi1_vec_sim = zeros(N,1);
    PHI_vec_sim = zeros(N,1);
    phi2_vec_sim = zeros(N,1);
    if 0
        disp_count = 0;
        disp_max = ceil(N/9.99);
        fprintf(' %% complete:  0%% ');
        for i = 1:N
            disp_count = disp_count+1;
            if disp_count == disp_max
                disp_count = 0;
                fprintf('%g%% ', round(i/N*100));
            end
            if i == N
                fprintf('100%% \n');
            end
            component_IX = MTR_IX_sim(i);
            [phi1_accept, PHI_accept, phi2_accept] = sample_orientation_from_ODF(ODF_best(component_IX), uniformODF);
            phi1_vec_sim(i) = phi1_accept;
            PHI_vec_sim(i) = PHI_accept;
            phi2_vec_sim(i) = phi2_accept;
        end
    else
        phi1_a = zeros(N,numel(ODF_best));
        PHI_a = zeros(N,numel(ODF_best));
        phi2_a = zeros(N,numel(ODF_best));
        for j = 1:numel(ODF_best)
            [phi1_accept, PHI_accept, phi2_accept] = sample_N_orientations_from_ODF(N, ODF_best(j), uniformODF);
            phi1_a(:,j) = phi1_accept;
            PHI_a(:,j) = PHI_accept;
            phi2_a(:,j) = phi2_accept;
        end
        for i = 1:N
            phi1_vec_sim(i,1) = phi1_a(i,MTR_IX_sim(i,1));
            PHI_vec_sim(i,1) = PHI_a(i,MTR_IX_sim(i,1));
            phi2_vec_sim(i,1) = phi2_a(i,MTR_IX_sim(i,1));
        end
    end
    fprintf('    complete.\n');
    
    fprintf('    writing to file...\n');
    save simulated_MTRs s_sim phi1_vec_sim PHI_vec_sim phi2_vec_sim
    clear phi1_accept PHI_accept phi2_accept phi1_a PHI_a phi2_a X_sim Y_sim
    fprintf('        complete.\n');
else
    load simulated_MTRs s_sim phi1_vec_sim PHI_vec_sim phi2_vec_sim 
end

%tmp = [s_sim, phi1_vec_sim, PHI_vec_sim, phi2_vec_sim];
%csvwrite('simulated_MTRs_csv.csv',tmp);

%%% plot IPF map

if 0
    fprintf('\ncalculating inverse pole figure...\n');
    
    figure(86)
    [X_fundamental,Y_fundamental,cmap1,cmap2,cmap3] = generate_inverse_pole_figure(phi1_vec_sim,PHI_vec_sim,phi2_vec_sim);
    clear X_fundamental Y_fundamental
    
    figure(87)
    cm = colormap('jet');
    cmin = 0;
    cmax = 1;
    caxis([cmin cmax])
    
    volume_element_side_length_sim = nx_sim * dx;
    fprintf('plotting...\n');
    a = zeros(ny,nx,3);
    k = 0;
    % zoi = 1;
    for zix = 1:nz
        for i = 1:nx   % columns
            for j = 1:ny   % rows
                k = k+1;
                if zix == zoi
                    ci1 = cmap1(k);
                    ci2 = cmap2(k);
                    ci3 = cmap3(k);
                    a(j, i,:) = [ci1, ci2, ci3];
                end
            end
        end
    end
    
    image([min(s_sim(:,1)) max(s_sim(:,1))], [min(s_sim(:,2)) max(s_sim(:,2))], a);
    set(gcf,'color','w');
    axis equal
    set(gca,'YDir','normal');
    imwrite(a,'sim_IPF_map.jpg','jpg');
    
    set(gcf,'color','w');
    hold on;
    plot(xlen*[0 1 1 0 0], ylen*[0 0 1 1 0],'k-');
    hold off;
    title('IPF Map');

else
    fprintf('\ncalculating inverse pole figure map...\n');
    a = view_IPF_map(s_sim, phi1_vec_sim, PHI_vec_sim, phi2_vec_sim);
    fprintf('    complete.\n');

    imwrite(a,'sim_IPF_map.jpg','jpg');
end

t = toc;
fprintf('total simulation compute time: %4.3f s\n',t);

%% try to plot 3D
% 
if 0
    [xx, yy, zz] = meshgrid(1:nx_sim, 1:ny_sim, 1:nz_sim);
    figure(3);
    
    hold on
    
    warp(xx, yy, ones(size(zz)), a)
    hold off
    set(gcf,'color','w');
    axis equal
end
