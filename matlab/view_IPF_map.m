% view_IPF_map.m

% Daniel M. Sparkman
% 07/03/2017
% Research

% view_IPF_map.m

% this code plots the IPF map

function [a] = view_IPF_map(s, phi1_vec, PHI_vec, phi2_vec, CI_mask)

%% setup

N = size(s,1);

x_vec = s(:,1);
y_vec = s(:,2);

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

nx = floor(xlen/dx + 1);
ny = floor(ylen/dy + 1);

if exist('CI_mask','var') == false
    CI_mask = false(ny,nx);
end

%% plot IPF map
%fprintf('\ncalculating inverse pole figure...\n');

figure(86);
cm = colormap('jet');
cmin = 0;
cmax = 1;
caxis([cmin cmax])

%fprintf('plotting...\n');
a = zeros(ny,nx,3);

[X,Y] = unit_triangle_IPF_coords(phi1_vec,PHI_vec,phi2_vec);
[cmap1,cmap2,cmap3] = IPF_colors(X,Y);

% if use_CI_mask
    cmap1(CI_mask) = cmap1(CI_mask)*0.0;
    cmap2(CI_mask) = cmap2(CI_mask)*0.0;
    cmap3(CI_mask) = cmap3(CI_mask)*0.0;
% end

nplot= min([1e4,N]);
ix = zeros(nplot,1);
rix = randperm(N);
len_CI_mask = numel(CI_mask);
kix = 0;
for aix = 1:nplot
    found_not_masked = false;
    while found_not_masked == false
        kix = kix + 1;
        if kix > N
            found_not_masked = true;
        else
            if rix(kix) <= len_CI_mask
                if CI_mask(rix(kix)) == false
                    ix(aix) = rix(kix);
                    found_not_masked = true;
                end
            end
        end
    end
end
        
scatter(X(ix),Y(ix),10,[cmap1(ix) cmap2(ix) cmap3(ix)],'fill' );

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

for k = 1:N
    ci1 = cmap1(k);
    ci2 = cmap2(k);
    ci3 = cmap3(k);
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
    a(j, i,:) = [ci1, ci2, ci3];
end

%% plot

figure(87);
image([min(x_vec) max(x_vec)], [min(y_vec) max(y_vec)], a);
set(gcf,'color','w');
axis equal
set(gca,'YDir','reverse');

set(gcf,'color','w');
hold on;
plot(xlen*[0 1 1 0 0], ylen*[0 0 1 1 0],'k-');
hold off;
title('IPF Map');
