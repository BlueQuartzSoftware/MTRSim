% generate_inverse_pole_figure.m

% Daniel M. Sparkman
% 02/14/2013
% Research
% generate_inverse_pole_figure.m

% function [X,Y] = generate_inverse_pole_figure(phi1,PHI,phi2)
function [X_fundamental,Y_fundamental,cmap1,cmap2,cmap3] = generate_inverse_pole_figure(phi1,PHI,phi2)

% hcp     fcc
crystal_system = 'hcp';


%%
% data = csvread('orientations_uniform.csv',2,0);
% phi1 = data(:,3);
% PHI = data(:,4);
% phi2 = data(:,5);

% data = csvread('orientations_textured.csv',2,0);
% phi1=data(:,3);
% PHI = data(:,4);
% phi2 = data(:,5);

% PHI = mod(PHI,pi/180*90);
% phi2 = mod(phi2,pi/180*60);

plot_on = true;

% set normal vector of specimen plane of interest in specimen coords
h = [0;
     0;
     1];

num_points = numel(phi1);
x = zeros(num_points,1);
y = zeros(num_points,1);
z = zeros(num_points,1);

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

    phi1s = symmetry_operators_euler_angles(:,1);
    PHIs = symmetry_operators_euler_angles(:,2);
    phi2s = symmetry_operators_euler_angles(:,3);
    
    theta_list = zeros(N,1);
    % r1_list = zeros(N,1);
    % r2_list = zeros(N,1);
    % r3_list = zeros(N,1);
    
    Xall = zeros(N*num_points,1);
    Yall = zeros(N*num_points,1);
    k = 0;
end

% x_fundamental = zeros(num_points,1);
% y_fundamental = zeros(num_points,1);
% z_fundamental = zeros(num_points,1);



%% step through each orientation

for i = 1:numel(phi1)
    c1 = cos(phi1(i));
    c2 = cos(phi2(i));
    C = cos(PHI(i));
    s1 = sin(phi1(i));
    s2 = sin(phi2(i));
    S = sin(PHI(i));
    
    % build rotation matrix
    g = [  c1*c2-s1*s2*C   s1*c2+c1*s2*C  s2*S ;
          -c1*s2-s1*c2*C  -s1*s2+c1*c2*C  c2*S ;
                s1*S           -c1*S        C  ];
    
    % rotate specimen plane normal vector into
    % crystal coord representation
    h_rot = g * h;
    
    % just examining half sphere because of symmetry over each side of the
    % plane
    if h_rot'*[0 0 -1]' < 0
        h_rot = -1 * h_rot;
    end
    
    % save unit sphere representation for this location's orientation
    x(i) = h_rot(1);
    y(i) = h_rot(2);
    z(i) = h_rot(3);
    
    if full_pole_figure
    % for this orientation, find symmetric orientation in fundamental zone
    x_symm = zeros(N,1);
    y_symm = zeros(N,1);
    z_symm = zeros(N,1);
    
    alpha_symm = zeros(N,1);
    beta_symm  = zeros(N,1);
    for j = 1:N
        c1s = cos(phi1s(j));
        c2s = cos(phi2s(j));
        Cs  = cos(PHIs(j));
        s1s = sin(phi1s(j));
        s2s = sin(phi2s(j));
        Ss  = sin(PHIs(j));
        
        % build crystal symmetry operator rotation matrix
        O_crystal = [  c1s*c2s-s1s*s2s*Cs   s1s*c2s+c1s*s2s*Cs  s2s*Ss ;
                      -c1s*s2s-s1s*c2s*Cs  -s1s*s2s+c1s*c2s*Cs  c2s*Ss ;
                              s1s*Ss               -c1s*Ss        Cs  ];
        
        % calculate new
%         g_rot_symm = g*O_crystal;   % I am not sure why this is wrong, but the equation below works!
        g_rot_symm = O_crystal *g;
        
        % rotate specimen plane normal vector into
        % crystal coord representation
        h_rot_symm = g_rot_symm * h;
        
        % just examining half sphere because of symmetry over each side of the
        % plane
        if h_rot_symm'*[0 0 -1]' < 0
            h_rot_symm = -1 * h_rot_symm;
        end
        
        k = k+1;
        Xall(k) = h_rot_symm(1) / (1-h_rot_symm(3));
        Yall(k) = h_rot_symm(2) / (1-h_rot_symm(3));
        
        alpha_symm(j) = acos(h_rot(3));
        beta_symm(j)  = atan2(h_rot(2),h_rot(1));

%         x_symm(j) = h_rot_symm(1);
%         y_symm(j) = h_rot_symm(2);
%         z_symm(j) = h_rot_symm(3);
    end
    
%     X_symm = x_symm./(1-z_symm);
%     Y_symm = y_symm./(1-z_symm);
%     
%     for j = 1:N
%         if X_symm(j) >= 0 && Y_symm(j) <= X_symm(j)*tan(pi/180*30) && Y_symm(j) >= 0
%             x_fundamental(i) = x_symm(j);
%             y_fundamental(i) = y_symm(j);
%             z_fundamental(i) = z_symm(j);
%         end
%     end
    end
	
end

% project unit sphere representation into stereographic plane (circle) 
X = x./(1-z);
Y = y./(1-z);

% % project unit sphere representation of fundamental zone into stereographic plane (circle) 
% X_fundamental = x_fundamental./(1-z_fundamental);
% Y_fundamental = y_fundamental./(1-z_fundamental);

X_fundamental = zeros(num_points,1);
Y_fundamental = zeros(num_points,1);
% only keep points in the standard stereographic triangle
switch crystal_system
    case 'hcp'
        k=0;
        for i = 1:num_points
            found_fundamental_flag = false;
            h = 0;
            while h < N && found_fundamental_flag == false
                h = h+1;
                j = (i-1)*N + h;
                if Xall(j) >= 0 && Yall(j) <= Xall(j)*tan(pi/180*30) && Yall(j) >= 0
                    k=k+1;
                    X_fundamental(k) = Xall(j);
                    Y_fundamental(k) = Yall(j);
                    found_fundamental_flag = true;
                end
                if h == N && found_fundamental_flag == false;
                    k = k+1;
                    fprintf('error in IPF calculation...for %gth orientation \n',i);
                end
            end
        end
    
    case 'fcc'
        
end

%% plot

if plot_on == true
%     figure(2)
    
    % colormap
    c = colormap('jet');
    cmin = 0;
    cmax = 1;
    caxis([cmin cmax])
    
    cmap1 = zeros(size(X_fundamental,1),1);
    cmap2 = zeros(size(X_fundamental,1),1);
    cmap3 = zeros(size(X_fundamental,1),1);
    set(gcf,'color','w');

    r = sqrt(X_fundamental.^2 +Y_fundamental.^2);
    beta = atan2(Y_fundamental, X_fundamental);
    % the color codings could be better, but probably good enough
    cmap1 = (1-r);
    cmap2 = r.*(1-beta/(pi/6));
    cmap3 = r.*(beta/(pi/6));
    max_c = max([cmap1 cmap2 cmap3],[],2);
    cmap1 = cmap1./max_c;
    cmap2 = cmap2./max_c;
    cmap3 = cmap3./max_c;
%     for i = 1:size(X_fundamental,1)
%         max_c = max([cmap1(i) cmap2(i) cmap3(i)]);
%         cmap1(i) = cmap1(i)/max_c;
%         cmap2(i) = cmap2(i)/max_c;
%         cmap3(i) = cmap3(i)/max_c;
%     end
%     figure(6)
%     hold on;
%     for i = 1:size(X_fundamental,1)
%         plot(X_fundamental(i),Y_fundamental(i),'o','MarkerSize',3,'MarkerFaceColor',[cmap1(i) cmap2(i) cmap3(i)],'MarkerEdgeColor',[cmap1(i) cmap2(i) cmap3(i)] );
%     end
%     colorbar
%     caxis([cmin cmax])
%     colorbar off
%     hold off;

    scatter(X_fundamental,Y_fundamental,10,[cmap1 cmap2 cmap3],'fill' );

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


%%
if plot_on == true
    if 0
%     figure(1)
    % plot original pole figure
    plot(X,Y,'bo','MarkerSize',3,'MarkerFaceColor','b')
    
    % generate outer circle
    theta = [0:0.1:2*pi,0];
    circle_x = cos(theta);
    circle_y = sin(theta);
    hold on;
    plot(circle_x,circle_y,'k-');
    hold off;
    set(gcf,'color','w');
    axis off;
    axis equal;
    text(-1.5, 1., '  original data', 'FontWeight', 'b', 'FontSize', 14);
    end
end


if plot_on == true && 0
%     figure(2)
    
    % colormap
    c = colormap('jet');
    cmin = 0;
    cmax = 1;
    caxis([cmin cmax])
%     
%     xc = [0 0.15 0.3 0.7 1 cos(pi/180*30)*[0.3 0.7 1] cos(pi/180*15)*[0.3 0.7 1]];
%     yc = [0 0.05 0   0   0 sin(pi/180*30)*[0.3 0.7 1] sin(pi/180*15)*[0.3 0.7 1]];
%     c1 = interp1(cmin:(cmax-cmin)/(64-1):cmax,c(:,1), [1 0.9 0.8 0.6 0.48 0.8 0.3 0 0.8 0.3 0.4]);
%     c2 = interp1(cmin:(cmax-cmin)/(64-1):cmax,c(:,2), [1 0.9 0.8 0.6 0.48 0.8 0.3 0 0.8 0.3 0.4]);
%     c3 = interp1(cmin:(cmax-cmin)/(64-1):cmax,c(:,3), [1 0.9 0.8 0.6 0.48 0.8 0.3 0 0.8 0.3 0.4]);

    xc = [0   0.25 0.45 0.6 1    cos(pi/180*30)*[0.3 0.5  0.7 0.9  1]  cos(pi/180*15)*[0.3 0.5 0.7 1]  ];
    yc = [0   0    0.05 0   0    sin(pi/180*30)*[0.3 0.4  0.7 0.9  1]  sin(pi/180*15)*[0.3 0.5 0.7 1]  ];
    c1 = [1   1    1    0.5 0.15                 1   0.75 0.5 0.25 0                   1.0 1.0 0   0.5 ];
    c2 = [0.0 0.5  1    1   0.5                  0.1 0.0  0   0.0  0                   0.5 1.0 1.0 0.5 ];
    c3 = [0.0 0    0    0.5 0.15                 0.2 0.75 1.0 1.0  1                   0.5 1.0 1.0 0.5 ];
    
    cmap1 = zeros(size(X_fundamental,1),1);
    cmap2 = zeros(size(X_fundamental,1),1);
    cmap3 = zeros(size(X_fundamental,1),1);
    set(gcf,'color','w');

    mu = 0.5;
    sigma = 0.5;
    theta = 0.2;
    
    rho = @(xi, xj, theta) exp(-(norm(xi-xj)/theta)^2);
    
    x_plot = [X_fundamental Y_fundamental];
    x = [xc' yc'];
    
    n = size(x,1);
    
    R = zeros(n,n);
    
    % correlation matrix
    for i = 1:n
        for j = 1:n
            R(i,j) = rho(x(i,:),x(j,:),theta);
        end
    end
    
    r_plot = zeros(n,1);
    for j = 1:size(x_plot,1)
        for i = 1:n
            r_plot(i) = rho(x(i,:),x_plot(j,:),theta);
        end
        w_plot = r_plot'*R^-1;
        cmap1(j) = mu + w_plot*(c1' - mu);
        cmap2(j) = mu + w_plot*(c2' - mu);
        cmap3(j) = mu + w_plot*(c3' - mu);
    end
    
    hold on;
    for i = 1:size(X_fundamental,1)

        bad_point_flag = false;
        
        % plot pole figure
        if cmap1(i) < 0
            cmap1(i) = 0;
            bad_point_flag = true;
        end
        if cmap1(i) > 1
            cmap1(i) = 1;
            bad_point_flag = true;
        end
        if cmap2(i) < 0
            cmap2(i) = 0;
            bad_point_flag = true;
        end
        if cmap2(i) > 1
            cmap2(i) = 1;
            bad_point_flag = true;
        end
        if cmap3(i) < 0
            cmap3(i) = 0;
            bad_point_flag = true;
        end
        if cmap3(i) > 1
            cmap3(i) = 1;
            bad_point_flag = true;
        end
            
        % show if a point had bad values or not...
        if bad_point_flag == true
            plot(X_fundamental(i),Y_fundamental(i),'o','MarkerSize',3,'MarkerFaceColor',[cmap1(i) cmap2(i) cmap3(i)],'MarkerEdgeColor',[cmap1(i) cmap2(i) cmap3(i)] );
        else
            plot(X_fundamental(i),Y_fundamental(i),'o','MarkerSize',3,'MarkerFaceColor',[cmap1(i) cmap2(i) cmap3(i)],'MarkerEdgeColor',[cmap1(i) cmap2(i) cmap3(i)] );
        end
    
    end
    colorbar
    caxis([cmin cmax])
    colorbar off

    hold off;

%     % plot pole figure
%     plot(X_fundamental,Y_fundamental,'o','MarkerSize',3,'MarkerFaceColor',[cmap1(i) cmap2(i))
%     
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
%     text(-.15, 0.6, '  fundamental', 'FontWeight', 'b', 'FontSize', 14);
    text(-.1, -.05, '  [0001]', 'FontWeight', 'b', 'FontSize', 14);
    text(0.9, -.05, '  [2-1-10]', 'FontWeight', 'b', 'FontSize', 14);
    text(0.8, 0.55, '  [10-10]', 'FontWeight', 'b', 'FontSize', 14);
end

if plot_on == true
    if 0
    if full_pole_figure
%         figure(3)
        % plot pole figure
        plot(Xall,Yall,'bo','MarkerSize',3,'MarkerFaceColor','b')
        
        % generate outer circle
        theta = [0:0.1:2*pi,0];
        circle_x = cos(theta);
        circle_y = sin(theta);
        hold on;
        plot(circle_x,circle_y,'k-');
        hold off;
        set(gcf,'color','w');
        axis off;
        axis equal;
        text(-1.05, 1.05, '  all orientation symmetries', 'FontWeight', 'b', 'FontSize', 14);
    end
    end
end

      