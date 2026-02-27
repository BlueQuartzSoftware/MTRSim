function [] = plot_PF(PF)

% plot_PF.m

% Daniel M. Sparkman
% 09/17/2013
% Research

% plot_PF.m

% this code plots the pole figure given the struct PF

%% set up

use_voronoi = false;
use_patch = true;

num_bins = PF.num_bins;

load uniformODF uniformODF
% temporarily put negligible but finite values in uniformODF to avoid
% divide by zero errors
for i = 1:length(uniformODF.ODFval(:,1))
    if uniformODF.ODFval(i,1) == 0
        uniformODF.ODFval(i,1) = 1e-15;
    end
end

% get a uniform PF (from uniform ODF if needed)
if 0	
    uniformPF = convert_ODF_to_PF(uniformODF.ODFval, uniformODF.ODFbins, uniformODF.phi1_bins, uniformODF.PHI_bins, uniformODF.phi2_bins);
    save uniformPF uniformPF
else
    load uniformPF uniformPF
end
PFval = PF.PFval;

% calc multiples of random (mrd) and put in 3rd column of PF struct
PFval(:,3) = PF.PFval(:,3) ./ uniformPF.PFval(:,3);

%% plot it

if use_voronoi
    % add zero outer ring for the actual outer ring to plot the voronoi cells
    r = 1+ 0.5/(num_bins/4);
    theta = 2*pi/(num_bins)/2 : 2*pi/(num_bins) : 2*pi - 2*pi/(num_bins)/2;
    xtmp = r*cos(theta');
    ytmp = r*sin(theta');
    PFval = [PFval; [xtmp ytmp zeros(num_bins,1)]];
end

if 0
    figure(1)
    hold on;
    for i = 1:size(PFval,1)
        if PFval(i,3) > 0
            plot3(PFval(i,1),PFval(i,2),PFval(i,3),'b.');
        end
    end
    hold off;
    xlabel('X');
    ylabel('Y');
    zlabel('I');
    set(gcf,'color','w');
else
    % plot it!
    x = PFval(:,1);
    y = PFval(:,2);
    
    if use_patch
        omit_center = true;
        % build patch
        c = {};
        v = [];
        % build the list of vertices
        beta = 0*pi/(num_bins-1) : 2*pi/(num_bins-1) : 2*pi;
        alpha = 0 : 2*pi/(num_bins-1) : pi/2;
        
        k = 0;
        l = 0;
        if omit_center
            % account for circle in center versus quads in outer rings
            alpha_list = [ alpha(2)*ones(1,numel(beta)) ];
            beta_list = beta;
            % convert to x,y for plotting
            r = tan(alpha_list / 2);
            x = r .* cos(beta_list);
            y = r .* sin(beta_list);
            % build vertex list for patch
            c_list = [];
            for ix = 1:numel(x)
                l = l+1;
                v(l,:) = [x(ix) y(ix)];
                c_list(ix) = l;
            end
            % build index list for patch
            k = k+1;
            c{k} = c_list;
        end
        for i = 2:numel(alpha)-1
            for j = 1:numel(beta)-1
                % account for triangles in center versus quads in outer rings
                if i > 1
                    alpha_list = [alpha(i) alpha(i+1) alpha(i+1) alpha(i)];
                    beta_list =  [beta(j) beta(j) beta(j+1) beta(j+1)];
                elseif omit_center == false
                    alpha_list = [alpha(i) alpha(i+1) alpha(i+1)];
                    beta_list = [beta(j) beta(j) beta(j+1)];
                end
                % convert to x,y for plotting
                r = tan(alpha_list / 2);
                x = r .* cos(beta_list);
                y = r .* sin(beta_list);
                % build vertex list for patch
                c_list = [];
                for ix = 1:numel(x)
                    l = l+1;
                    v(l,:) = [x(ix) y(ix)];
                    c_list(ix) = l;
                end
                % build index list for patch
                k = k+1;
                c{k} = c_list;
            end
        end
    end
    
    if use_voronoi
        %     figure(1)
        [TRI] = delaunay(x,y);
        %     voronoi(x,y)
        [v,c] = voronoin([x y]);
    end
    
    cm = colormap('jet');
    cmin = 0;
    if omit_center
        dist_sq = sum(([PFval(:,1), PFval(:,2)]).^2,2);
        center_ix = sqrt(dist_sq) < tan(alpha(2) / 2);
        center_PFval = mean(PFval(center_ix,3));
        cmax = max([center_PFval; PFval(~center_ix,3)]); %cmax = 3.86;% !!!
    else
        cmax = max(PFval(:,3));
    end
    
    caxis([cmin cmax])
    opengl software

    
    % plot the pole figure intensity values according to colorbar
    for i = 1:length(c)
        if omit_center 
            if i == 1
                ci1 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,1),center_PFval,'linear');
                ci2 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,2),center_PFval,'linear');
                ci3 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,3),center_PFval,'linear');
                if center_PFval > cmax
                    ci1 = 1;
                    ci2 = 1;
                    ci3 = 1;
                end
            else
                ix = i + numel(beta)-1;
                ci1 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,1),PFval(ix,3),'linear');
                ci2 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,2),PFval(ix,3),'linear');
                ci3 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,3),PFval(ix,3),'linear');
                if PFval(ix,3) > cmax
                    ci1 = 1;
                    ci2 = 1;
                    ci3 = 1;
                end
            end
        else
                ci1 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,1),PFval(i,3),'linear');
                ci2 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,2),PFval(i,3),'linear');
                ci3 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,3),PFval(i,3),'linear');
        end
        
        tcolor = [ci1 ci2 ci3];  % pole figure intensity
        
        if use_voronoi 
            patch( v(c{i},1) , v(c{i},2) , tcolor, 'EdgeAlpha',0 );
        elseif use_patch
            if i > 1
                % approx centroid of patch to color
                x_avg = mean( v(c{i},1) );
                y_avg = mean( v(c{i},2) );
                % find index of PFval that is closet to the patch centroid
                dist_sq = sum(([PFval(:,1) - x_avg, PFval(:,2) - y_avg]).^2,2);
                [val,ix] = min(dist_sq);
                % determine color of patch
                ci1 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,1),PFval(ix,3),'linear');
                ci2 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,2),PFval(ix,3),'linear');
                ci3 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,3),PFval(ix,3),'linear');
                if PFval(ix,3) > cmax
                    ci1 = 1;
                    ci2 = 1;
                    ci3 = 1;
                end
                tcolor = [ci1 ci2 ci3];  % pole figure intensity
            else
                ci1 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,1),center_PFval,'linear');
                ci2 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,2),center_PFval,'linear');
                ci3 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,3),center_PFval,'linear');
                tcolor = [ci1 ci2 ci3];  % pole figure intensity
            end
            % plot colored patch
            patch( v(c{i},1) , v(c{i},2) , tcolor, 'EdgeAlpha',0 );
        else
            hold on;
            plot(x(i),y(i),'o','MarkerFaceColor',tcolor,'MarkerEdgeColor',tcolor);
            hold off;
        end
    end
    
    colorbar
    caxis([cmin cmax])
    colorbar off
    
    set(gcf,'color','w');
    axis equal;
    axis off;
end