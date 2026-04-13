function [x] = norminv_code(p)

% norminv_code.m

% Daniel M. Sparkman
% 08/28/2013
% Research

% norminv_code.m

% this code computes the inverse normal cdf evaluated for p
% Uses erfcinv (built-in) so the Statistics Toolbox is not required.
z = -sqrt(2) * erfcinv( 2*p );
mu = 0;
sig = 1;
x = mu + sig * z;