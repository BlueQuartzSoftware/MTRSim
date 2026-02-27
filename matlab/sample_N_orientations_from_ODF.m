function [phi1_accept, PHI_accept, phi2_accept] = sample_N_orientations_from_ODF(N, ODF, uniformODF)

% sample_orientation_from_ODF.m

% Daniel M. Sparkman
% 08/27/2013
% Research

% sample_orientation_from_ODF.m

degree_spacing = 180/pi*(ODF.phi1_bins(2)-ODF.phi1_bins(1));
bin_scaling_factor = degree_spacing/180*pi;

u0 = rand(1,N);

% construct the bins
edges = [0, cumsum(ODF.ODFval(:,1))'];
s = edges(end);
if abs(s - 1) > eps
    edges = edges * (1 / s);
end

% draw bins
c = histc(u0, edges);
ce = c(end);
c = c(1:end-1);
c(end) = c(end) + ce;

% extract samples

xv = find(c);

if numel(xv) == N  % each value is sampled at most once
    ix = xv;
else                % some values are sampled more than once
    xc = c(xv);
    d = zeros(1, N);
    dv = [xv(1), diff(xv)];
    dp = [1, 1 + cumsum(xc(1:end-1))];
    d(dp) = dv;
    ix = cumsum(d);
end

% randomly permute the sample's order
ix = ix(randperm(N));

% ix = find(cumsum(ODF.ODFval(:,1)) >= u0,1,'first');

phi1_accept = uniformODF.ODFbins(ix,1) + bin_scaling_factor*(rand(N,1)-0.5);
PHI_accept  = uniformODF.ODFbins(ix,2) + bin_scaling_factor*(rand(N,1)-0.5);
phi2_accept = uniformODF.ODFbins(ix,3) + bin_scaling_factor*(rand(N,1)-0.5);

