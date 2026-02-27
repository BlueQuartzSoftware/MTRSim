function [] = plot_AR(min_thresholds, max_thresholds, custom_map)
% plot_AR.m

% Daniel M. Sparkman
% 09/18/2013
% Research

% plot_AR.m

% this code takes the thresholds in an assignment rule and plots them
% block assignments are assumed

num_components = size(min_thresholds, 1);
num_gaussians  = size(min_thresholds, 2);

switch num_gaussians
    case 1
% % % %         % build plot location list
% % % %         num_pts = num_components;
% % % %         s = zeros(num_pts,num_gaussians);
% % % %         for k = 1:num_pts
% % % %             for j = 1:num_gaussians
% % % %                 s(k,j) = (min_thresholds(k,j) + max_thresholds(k,j))/2 ;
% % % %             end
% % % %         end
            
        % assign color to locations
        AR_IX = 1:num_components;
        
        % plot
        cm = colormap(custom_map);
        cmin = min(AR_IX);
        cmax = max(AR_IX);
        caxis([cmin cmax])
        hold on;
        for i = 1:num_components
            ci1 = cm(AR_IX(i),1);
            ci2 = cm(AR_IX(i),2);
            ci3 = cm(AR_IX(i),3);
%             ci1 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,1),AR_IX(i),'linear');
%             ci2 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,2),AR_IX(i),'linear');
%             ci3 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,3),AR_IX(i),'linear');
            
            tcolor(1,i,1:3) = [ci1 ci2 ci3];  % microtexture assignment
            
            v = [min_thresholds(i,1) 0;
                 max_thresholds(i,1) 0;
                 max_thresholds(i,1) 1;
                 min_thresholds(i,1) 1];
            
            patch( v(:,1) , v(:,2) ,tcolor(1,i,1:3), 'EdgeAlpha',0);
        end
        hold off;
        axis equal
        colorbar
        caxis([cmin cmax])
        set(gcf,'color','w');
        axis off;
        title('Assignment Rule');
        
    case 2
% % % %         % build plot location list
% % % %         num_pts = num_components;
% % % %         s = zeros(num_pts,num_gaussians);
% % % %         for k = 1:num_pts
% % % %             for j = 1:num_gaussians
% % % %                 s(k,j) = (min_thresholds(k,j) + max_thresholds(k,j))/2 ;
% % % %             end
% % % %         end
            
        % assign color to locations
        AR_IX = 1:num_components;
        
        % plot
        cm = colormap(custom_map);
        cmin = min(AR_IX);
        cmax = max(AR_IX);
        caxis([cmin cmax])
        opengl software

%         hold on;
        for i = 1:num_components
            ci1 = cm(AR_IX(i),1);
            ci2 = cm(AR_IX(i),2);
            ci3 = cm(AR_IX(i),3);
%             ci1 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,1),AR_IX(i),'linear');
%             ci2 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,2),AR_IX(i),'linear');
%             ci3 = interp1(cmin:(cmax-cmin)/(64-1):cmax,cm(:,3),AR_IX(i),'linear');
            
            tcolor(1,i,1:3) = [ci1 ci2 ci3];  % microtexture assignment
            
            v = [min_thresholds(i,1) min_thresholds(i,2);
                 max_thresholds(i,1) min_thresholds(i,2);
                 max_thresholds(i,1) max_thresholds(i,2);
                 min_thresholds(i,1) max_thresholds(i,2)];
            
            patch( v(:,1) , v(:,2) ,tcolor(1,i,1:3), 'EdgeAlpha',0);
        end
%         hold off;
        colorbar
        caxis([cmin cmax])
                
        set(gcf,'color','w');
        axis equal;
        axis off;
        title('Assignment Rule');
%         set(gcf, 'renderer', 'zbuffer');
    case 4
        fprintf(' unable to plot assignment rule for 4 gaussians...\n');
end
        