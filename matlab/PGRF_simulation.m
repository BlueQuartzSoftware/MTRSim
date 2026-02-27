function [MTR_IX, Z_all] = PGRF_simulation(P, theta_list, nugvar, dx, dy, dz, nx, ny, nz)
% PGRF_simulation.m
%--------------------------------------------------------------------------
% Daniel M. Sparkman
% 08/28/2013
% Research
% 
%--------------------------------------------------------------------------

% PGRF_simulation.m

% this code simulates assignments using the plurigaussian simulation technique

%% set up

% anisotropic     isotropic    radial
correlation_function_selected = 'anisotropic';

%    exp    gaussian    spherical    Matern
% !!!!! the Matern correlation function has not been implemented yet!!!!
% !! Spherical is not an option because it is not separable!!
corr_func_name = 'exp'; %'gaussian';

% periodic        nonperiodic
boundary_conditions = 'nonperiodic';

% spot       stationary
mean_function_selected = 'stationary';

N = nx*ny*nz; % number of points in underlying field to simulate
volume_element_side_length = max([nx*dx, ny*dy, nz*dz]);
[num_gaussians, AR_index, min_thresholds, max_thresholds, P_min_thresholds, P_max_thresholds] = select_AR(P);
mu_list = zeros(num_gaussians,1);

% mean of the random field for transformation purposes (usually just set to zero)
mu_field = 0;
% standard deviation of the random field for transformations purposes (usually just set to one)
sig_ensemble = 1;

switch boundary_conditions
    case 'nonperiodic'
        switch corr_func_name
            case 'exp'
                 rho = @(x,theta) exp( -(abs(x)/theta) );
            case 'gaussian'
                 rho = @(x,theta) exp( -(x/theta).^2 );
        end
    case 'periodic'
        rho = @(x,theta) exp( -(  min([abs(x), abs(x-volume_element_side_length), abs(x+volume_element_side_length)  ])/theta).^2 );
end

switch correlation_function_selected
    case 'isotropic'
        theta_list = [theta_list(:,1) theta_list(:,1) theta_list(:,1)];
end

fprintf('mean function:          %s\n', mean_function_selected);
fprintf('correlation function:   %s\n', correlation_function_selected);
fprintf('correlation model:      %s\n', corr_func_name);
fprintf('boundary conditions:    %s\n', boundary_conditions);
fprintf('\n');


%% build parameters of PGRF

Z_all = zeros(N,num_gaussians);


%% simulate underlying gaussians

for h = 1:num_gaussians
    fprintf('\ngaussian Y%g\n',h);
    
        % calc mean of realization of microtexture probabilities from P_A_mean
        mu_const = mu_list(h);
        
        % correlation parameters
        theta = theta_list(h,:);
        
        fprintf('building mean vector...\n');
        switch mean_function_selected
            case 'spot'
                
                % set location of the center of the spot
                spot_center = [0.5, 0.5, 0.5];
                
                % build mean vector
                mu_vec = zeros(N,1);
                for i = 1:N
                    if norm( [1 3].*([x(i) y(i) z(i)] - spot_center) ) < 0.3
                        mu_vec(i) = 1 - mu_const;
                    else
                        mu_vec(i) = mu_const;
                    end
                end
                % set standard deviation of random field
                sig_ensemble = 1;
                
            case 'stationary'
                mu_vec = mu_const*ones(N,1);
                sig_ensemble = 1;
        end

        fprintf('simulating underlying Gaussian random field...\n');
        
        if 1 
            tic; Z = gp_generator(dx, dy, dz, rho, theta, nx, ny, nz); toc
        else
            Z = 1;
        end
        Z_all(:,h) = Z(:)+mu_vec;

end

%% apply assignment rule (rock type rule)


fprintf('\napplying assignment rule...\n');

fprintf('data: \n');
num_components = numel(P);   % facies

for j = 1:num_components
    fprintf('P%g  =  %3.2f\n',j,P(j));
end

% figure(105);
% custom_map = [0.1 0.1 0.1;
%     1 0.9 0;
%     1 0 1;
%     0 1 0.1;
%     0 0 1];
% plot_AR(P_min_thresholds, P_max_thresholds, custom_map);

tic; MTR_IX = eval_AR(Z_all, min_thresholds, max_thresholds); toc

fprintf('simulated:\n');
e_P = zeros(num_components,1);
for j = 1:num_components
    e_P(j) = sum(MTR_IX==j)/N;
    fprintf('P%g  =  %3.2f\n',j,e_P(j));
end


