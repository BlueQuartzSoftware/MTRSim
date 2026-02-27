function [Assignments] = eval_AR(Z, min_thresholds, max_thresholds)

% eval_AR.m

% Daniel M. Sparkman
% 9/19/2013
% Research
% eval_AR.m

% this code evaluates an Assignment Rule given a set of plurigaussian
% simulations and determines the corresponding Assignments

%% set up

N = size(Z,1);
num_gaussians = size(Z,2);
num_components = size(min_thresholds, 1);

Assignments = zeros(N,1);
for i = 1:N
    component_list = 1:num_components;
    for j = 1:num_components
        for k = 1:num_gaussians
            if (Z(i,k) < min_thresholds(j,k)) || (Z(i,k) > max_thresholds(j,k))
                component_list(j) = 0;
            end
        end
    end
    Assignments(i) = max(component_list);
end
