function [phi1_accept, PHI_accept, phi2_accept] = sample_orientation_from_ODF(ODF, uniformODF)

% sample_orientation_from_ODF.m

% Daniel M. Sparkman
% 08/27/2013
% Research
% sample_orientation_from_ODF.m

use_acceptance_rejection = false;
degree_spacing = 180/pi*(ODF.phi1_bins(2)-ODF.phi1_bins(1));

if use_acceptance_rejection 
% this code uses acceptance rejection sampling to sample from a target ODF
% using a uniformly random texture ODF (Shoemake algorithm)

% determine scaling factor for acceptance probability calculation
M = max( ODF.ODFval(:,1) ./ uniformODF.ODFval(:,1) );

% number of samples to evaluate at a time
k = 1e2;
f_g = zeros(k,1);
g_g = zeros(k,1);

accept_sample = false;
sample_count = 0;
jumping_width = 0;
while accept_sample == false
    u0 = rand(k,1);
    
    % generate quaternion
    u1 = rand(k,1);
    u2 = rand(k,1);
    u3 = rand(k,1);
    q0 = sqrt(1-u1).*sin(2*pi*u2);
    q1 = sqrt(1-u1).*cos(2*pi*u2);
    q2 = sqrt(u1).*sin(2*pi*u3);
    q3 = sqrt(u1).*cos(2*pi*u3);
    
    % build rotation matrix
    A11 = 1-2*(q2.*q2 + q3.*q3);
    A12 = 2*q1.*q2 - 2*q0.*q3;
    A13 = 2*q1.*q3 + 2*q0.*q2;
%     A21 = 2*q1.*q2+2*q0.*q3;
%     A22 = 1-2*(q1.*q1 + q3.*q3);
    A23 = 2*q2.*q3 - 2*q0.*q1;
    A31 = 2*q1.*q3-2*q0.*q2;
    A32 = 2*q2.*q3 + 2*q0.*q1;
    A33 = 1-2*(q1.*q1 + q2.*q2);
    
    % R = [A11 A12 A13;
    %     A21 A22 A23;
    %     A31 A32 A33]'
    
    PHI = acos(A33);
    if PHI > 0
        phi1 = atan2(A31,(-A32));
        phi2 = atan2(A13,A23);
    elseif PHI == 0
        phi1 = atan2(A12,A11);
        phi2 = 0;
    else
        phi1 = atan2(A31,(-A32));
        phi2 = atan2(A13,A23);
    end
    
    % shift negative values of the euler angles to [0,2*pi]
    phi1(phi1 < 0) = phi1(phi1 < 0)+2*pi;
    phi2(phi2 < 0) = phi2(phi2 < 0)+2*pi;
    
    for i = 1:k
    % evaluate target distribution pdf (ODF)
    f_g(i) = ODF_discrete_eval(phi1(i),PHI(i),phi2(i),ODF.ODFval);
    
    % evaluate sampling distribution pdf (uniform)
    g_g(i) = ODF_discrete_eval(phi1(i),PHI(i),phi2(i),uniformODF.ODFval);
    end
    P_B = f_g ./ M ./ g_g;
    
    % accept or reject sample
    if any(u0 < P_B)
        sample_count = sample_count+1;
        if sample_count > jumping_width
            sample_count = 0;
            accept_sample = true;
        end
    else
        % reject sample and go through loop again
    end
    
end

ix = find(u0 < P_B,1,'first');
phi1_accept = phi1(ix);
PHI_accept = PHI(ix);
phi2_accept = phi2(ix);

else
    u0 = rand();
    bin_scaling_factor = degree_spacing/180*pi;
    
    ix = find(cumsum(ODF.ODFval(:,1)) >= u0,1,'first');
    
    phi1_accept = uniformODF.ODFbins(ix,1) + bin_scaling_factor*(rand()-0.5);
    PHI_accept  = uniformODF.ODFbins(ix,2) + bin_scaling_factor*(rand()-0.5);
    phi2_accept = uniformODF.ODFbins(ix,3) + bin_scaling_factor*(rand()-0.5);
end
