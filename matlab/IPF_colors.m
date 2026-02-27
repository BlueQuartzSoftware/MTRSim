% IPF_colors.m

% Daniel M. Sparkman
% 07/05/2017
% Research

% IPF_colors.m

% this code takes as inputs the x- and y- coordinates in the standard 
% stereoscopic unit triangle and outputs the colormap for an IPF

function [cmap1,cmap2,cmap3] = IPF_colors(X,Y)

%% setup

N = size(X,1);

plot_on = false;

%% plot

%     figure(2)
    
    % colormap
    c = colormap('jet');
    cmin = 0;
    cmax = 1;
    caxis([cmin cmax])
    
    r = sqrt(X.^2 +Y.^2);
    beta = atan2(Y, X);
    % the color codings could be better, but probably good enough
    cmap1 = (1-r);
    cmap2 = r.*(1-beta/(pi/6));
    cmap3 = r.*(beta/(pi/6));
    max_c = max([cmap1 cmap2 cmap3],[],2);
    cmap1 = cmap1./max_c;
    cmap2 = cmap2./max_c;
    cmap3 = cmap3./max_c;

if plot_on == true
    set(gcf,'color','w');
    scatter(X,Y,10,[cmap1 cmap2 cmap3],'fill' );

    % generate outer circle
    theta = [0:0.1:pi/6,pi/6];
    circle_x = cos(theta);
    circle_y = sin(theta);
    hold on;
    plot(circle_x,circle_y,'k-');
    plot([0 1],[0 0],'k-');
    plot([0 cos(pi/6)],[0 sin(pi/6)],'k-');
    hold off;
    set(gcf,'color','w');
    axis off;
    axis equal;
    text(-.1, -.05, '  [0001]', 'FontWeight', 'b', 'FontSize', 14);
    text(0.9, -.05, '  [2-1-10]', 'FontWeight', 'b', 'FontSize', 14);
    text(0.8, 0.55, '  [10-10]', 'FontWeight', 'b', 'FontSize', 14);
end


