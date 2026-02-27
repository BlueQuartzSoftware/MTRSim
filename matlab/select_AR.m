function [num_gaussians, AR_index, min_thresholds, max_thresholds, P_min_thresholds, P_max_thresholds] = select_AR(P)

% select_AR.m

% Daniel M. Sparkman
% 08/27/2013
% Research

% select_AR.m

num_components = length(P);
num_gaussians = num_components -1;
AR_index = num_components;
P_tol = 1e-6;

switch AR_index
    case 1    % 3 categories all touching
        P_min_thresholds(1,:) = [  0  0 ];
        P_max_thresholds(1,:) = [ P(1)  1 ];
        P_min_thresholds(2,:) = [ P(1)  0 ];
        P_max_thresholds(2,:) = [  1  P(2)/(1-P(1)) ];
        P_min_thresholds(3,:) = [ P(1)  P(2)/(1-P(1)) ];
        P_max_thresholds(3,:) = [  1  1 ];
        
        % determine thresholds for first component
        s = [-inf inf];
        t = [-inf inf];
        q = 0;
        qmin = -6;
        qmax = 6;
        Poi = P(1);
        P_tmp = qsimvn(20,[1 0; 0 1],[s(1) t(1)],[q t(2)]);
        while abs(Poi - P_tmp) > P_tol
            P_tmp = qsimvn(20,[1 0; 0 1],[s(1) t(1)],[q t(2)]);
            if P_tmp < Poi
                qmin = q;
                q = (qmax + qmin)/2;                
            elseif P_tmp > Poi
                qmax = q;
                q = (qmax + qmin)/2;
            end
        end
        s(2) = q;
        min_thresholds(1,:) = [ s(1) t(1) ];
        max_thresholds(1,:) = [ s(2) t(2) ];
        
        % determine thresholds for second component
        s = [max_thresholds(1,1) inf];
        t = [-inf inf];
        q = 0;
        qmin = -6;
        qmax = 6;
        Poi = P(2);
        P_tmp = qsimvn(20,[1 0; 0 1],[s(1) t(1)],[t(1) t(2)]);
        while abs(Poi - P_tmp) > P_tol
            P_tmp = qsimvn(20,[1 0; 0 1],[s(1) t(1)],[s(2) q]);
            if P_tmp < Poi
                qmin = q;
                q = (qmax + qmin)/2;                
            elseif P_tmp > Poi
                qmax = q;
                q = (qmax + qmin)/2;
            end
        end
        t(2) = q;
        min_thresholds(2,:) = [ s(1) t(1) ];
        max_thresholds(2,:) = [ s(2) t(2) ];
        
        % determine thresholds for third component
        s = [max_thresholds(1,1) inf];
        t = [max_thresholds(2,2) inf];
        min_thresholds(3,:) = [ s(1) t(1) ];
        max_thresholds(3,:) = [ s(2) t(2) ];
        
    case 2    % 5 categories but not all touching
        % num rows = num_components
        % num_cols = num_gaussians
        P_min_thresholds(1,:) = [  0  0  ];
        P_max_thresholds(1,:) = [ P(1)+P(2)  P(1)/(P(1)+P(2))  ];
        P_min_thresholds(2,:) = [  0  P(1)/(P(1)+P(2)) ];
        P_max_thresholds(2,:) = [  P(1)+P(2)  1 ];
        P_min_thresholds(3,:) = [ P(1)+P(2)  0 ];
        P_max_thresholds(3,:) = [  1  P(3)/(P(3)+P(4)+P(5)) ];
        P_min_thresholds(4,:) = [ P(1)+P(2)  P(3)/(P(3)+P(4)+P(5)) ];
        P_max_thresholds(4,:) = [  1  (P(3)+P(4))/(P(3)+P(4)+P(5)) ];
        P_min_thresholds(5,:) = [ P(1)+P(2)  (P(3)+P(4))/(P(3)+P(4)+P(5)) ];
        P_max_thresholds(5,:) = [  1  1 ];
        
        % determine thresholds for first component
        s = [-inf inf];
        t = [-inf inf];
        % find z1 first
        q = 0;
        qmin = -6;
        qmax = 6;
        Poi = P(1)+P(2);
        P_tmp = qsimvn(20,[1 0; 0 1],[s(1) t(1)],[q t(2)]);
        while abs(Poi - P_tmp) > P_tol
            P_tmp = qsimvn(20,[1 0; 0 1],[s(1) t(1)],[q t(2)]);
            if P_tmp < Poi
                qmin = q;
                q = (qmax + qmin)/2;                
            elseif P_tmp > Poi
                qmax = q;
                q = (qmax + qmin)/2;
            end
        end
        s(2) = q;
        % find z2 next
        q = 0;
        qmin = -6;
        qmax = 6;
        Poi = P(1);
        P_tmp = qsimvn(20,[1 0; 0 1],[s(1) t(1)],[s(2) t(2)]);
        while abs(Poi - P_tmp) > P_tol
            P_tmp = qsimvn(20,[1 0; 0 1],[s(1) t(1)],[s(2) q]);
            if P_tmp < Poi
                qmin = q;
                q = (qmax + qmin)/2;                
            elseif P_tmp > Poi
                qmax = q;
                q = (qmax + qmin)/2;
            end
        end
        t(2) = q;
        min_thresholds(1,:) = [ s(1) t(1) ];
        max_thresholds(1,:) = [ s(2) t(2) ];
        
        % determine thresholds for second component
        s = [-inf max_thresholds(1,1)];
        t = [max_thresholds(1,2) inf];
        min_thresholds(2,:) = [ s(1) t(1) ];
        max_thresholds(2,:) = [ s(2) t(2) ];
        
        % determine thresholds for third component
        s = [max_thresholds(1,1) inf];
        t = [-inf inf];
        q = 0;
        qmin = -6;
        qmax = 6;
        Poi = P(3);
        P_tmp = qsimvn(20,[1 0; 0 1],[s(1) t(1)],[s(1) t(2)]);
        while abs(Poi - P_tmp) > P_tol
            P_tmp = qsimvn(20,[1 0; 0 1],[s(1) t(1)],[s(2) q]);
            if P_tmp < Poi
                qmin = q;
                q = (qmax + qmin)/2;                
            elseif P_tmp > Poi
                qmax = q;
                q = (qmax + qmin)/2;
            end
        end
        t(2) = q;
        min_thresholds(3,:) = [ s(1) t(1) ];
        max_thresholds(3,:) = [ s(2) t(2) ];
        
        % determine thresholds for fourth component
        s = [min_thresholds(3,1) inf];
        t = [max_thresholds(3,2) inf];
        q = 0;
        qmin = -6;
        qmax = 6;
        Poi = P(4);
        P_tmp = qsimvn(20,[1 0; 0 1],[s(1) t(1)],[t(1) t(2)]);
        while abs(Poi - P_tmp) > P_tol
            P_tmp = qsimvn(20,[1 0; 0 1],[s(1) t(1)],[s(2) q]);
            if P_tmp < Poi
                qmin = q;
                q = (qmax + qmin)/2;                
            elseif P_tmp > Poi
                qmax = q;
                q = (qmax + qmin)/2;
            end
        end
        t(2) = q;
        min_thresholds(4,:) = [ s(1) t(1) ];
        max_thresholds(4,:) = [ s(2) t(2) ];
        
        % determine thresholds for fifth component
        s = [min_thresholds(4,1) inf];
        t = [max_thresholds(4,2) inf];
        min_thresholds(5,:) = [ s(1) t(1) ];
        max_thresholds(5,:) = [ s(2) t(2) ];
        
    case 5    % 5 categories all touching
        num_gaussians = num_components - 1;
        % P_min_thresholds and P_max_thresholds only used for plotting...
        % num rows = num_components
        % num_cols = num_gaussians
        P_min_thresholds(1,:) = [  0  0  0  0 ];
        P_max_thresholds(1,:) = [  1  1  1  1 ];
        P_min_thresholds(2,:) = [  0  0  0  0 ];
        P_max_thresholds(2,:) = [  1  1  1  1 ];
        P_min_thresholds(3,:) = [  0  0  0  0 ];
        P_max_thresholds(3,:) = [  1  1  1  1 ];
        P_min_thresholds(4,:) = [  0  0  0  0 ];
        P_max_thresholds(4,:) = [  1  1  1  1 ];
        P_min_thresholds(5,:) = [  0  0  0  0 ];
        P_max_thresholds(5,:) = [  1  1  1  1 ];
        
        % determine thresholds for first component
        s = -inf * ones(1, num_gaussians);
        t =  inf * ones(1, num_gaussians);
        % find z1 first
        q = 0;
        qmin = -6;
        qmax = 6;
        Poi = P(1);
        goi = 1;
        s_tmp = s;
        t_tmp = t; t_tmp(goi) = q;
        P_tmp = qsimvn(20,eye(num_gaussians),s_tmp,t_tmp);
        while abs(Poi - P_tmp) > P_tol
            t_tmp(goi) = q;
            P_tmp = qsimvn(20,eye(num_gaussians),s_tmp,t_tmp);
            if P_tmp < Poi
                qmin = q;
                q = (qmax + qmin)/2;                
            elseif P_tmp > Poi
                qmax = q;
                q = (qmax + qmin)/2;
            end
        end
        t(goi) = q;
        min_thresholds(1,:) = s;
        max_thresholds(1,:) = t;
        
        % determine thresholds for second component
        coi = 2;
        goi = 2;
        s = -inf * ones(1, num_gaussians);
        for i = 1:coi-1
            s(i) = max_thresholds(i,i);
        end
        t =  inf * ones(1, num_gaussians);
        q = 0;
        qmin = -6;
        qmax = 6;
        Poi = P(coi);
        s_tmp = s;
        t_tmp = t; t_tmp(goi) = q;
        P_tmp = qsimvn(20,eye(num_gaussians),s_tmp,t_tmp);
        while abs(Poi - P_tmp) > P_tol
            t_tmp(goi) = q;
            P_tmp = qsimvn(20,eye(num_gaussians),s_tmp,t_tmp);
            if P_tmp < Poi
                qmin = q;
                q = (qmax + qmin)/2;                
            elseif P_tmp > Poi
                qmax = q;
                q = (qmax + qmin)/2;
            end
        end
        t(goi) = q;
        min_thresholds(coi,:) = s;
        max_thresholds(coi,:) = t;
        
        % determine thresholds for third component
        coi = 3;
        goi = 3;
        s = -inf * ones(1, num_gaussians);
        for i = 1:coi-1
            s(i) = max_thresholds(i,i);
        end
        t =  inf * ones(1, num_gaussians);
        q = 0;
        qmin = -6;
        qmax = 6;
        Poi = P(coi);
        s_tmp = s;
        t_tmp = t; t_tmp(goi) = q;
        P_tmp = qsimvn(20,eye(num_gaussians),s_tmp,t_tmp);
        while abs(Poi - P_tmp) > P_tol
            t_tmp(goi) = q;
            P_tmp = qsimvn(20,eye(num_gaussians),s_tmp,t_tmp);
            if P_tmp < Poi
                qmin = q;
                q = (qmax + qmin)/2;                
            elseif P_tmp > Poi
                qmax = q;
                q = (qmax + qmin)/2;
            end
        end
        t(goi) = q;
        min_thresholds(coi,:) = s;
        max_thresholds(coi,:) = t;
        
        % determine thresholds for fourth component
        coi = 4;
        goi = 4;
        s = -inf * ones(1, num_gaussians);
        for i = 1:coi-1
            s(i) = max_thresholds(i,i);
        end
        t =  inf * ones(1, num_gaussians);
        q = 0;
        qmin = -6;
        qmax = 6;
        Poi = P(coi);
        s_tmp = s;
        t_tmp = t; t_tmp(goi) = q;
        P_tmp = qsimvn(20,eye(num_gaussians),s_tmp,t_tmp);
        while abs(Poi - P_tmp) > P_tol
            t_tmp(goi) = q;
            P_tmp = qsimvn(20,eye(num_gaussians),s_tmp,t_tmp);
            if P_tmp < Poi
                qmin = q;
                q = (qmax + qmin)/2;                
            elseif P_tmp > Poi
                qmax = q;
                q = (qmax + qmin)/2;
            end
        end
        t(goi) = q;
        min_thresholds(coi,:) = s;
        max_thresholds(coi,:) = t;
        
        % determine thresholds for fifth component
        coi = 5;
        goi = 5;
        s = -inf * ones(1, num_gaussians);
        for i = 1:coi-1
            s(i) = max_thresholds(i,i);
        end
        t =  inf * ones(1, num_gaussians);
        min_thresholds(5,:) = s;
        max_thresholds(5,:) = t;
        
    otherwise    % n categories all touching
        num_gaussians = num_components - 1;
        % P_min_thresholds and P_max_thresholds only used for plotting...
        % num rows = num_components
        % num_cols = num_gaussians
        P_min_thresholds = zeros(num_components,num_gaussians);
        P_max_thresholds = zeros(num_components,num_gaussians);
               
        % determine thresholds for first component
        s = -inf * ones(1, num_gaussians);
        t =  inf * ones(1, num_gaussians);
        % find z1 first
        q = 0;
        qmin = -6;
        qmax = 6;
        Poi = P(1);
        goi = 1;
        s_tmp = s;
        t_tmp = t; t_tmp(goi) = q;
        P_tmp = qsimvn(20,eye(num_gaussians),s_tmp,t_tmp);
        while abs(Poi - P_tmp) > P_tol
            t_tmp(goi) = q;
            P_tmp = qsimvn(20,eye(num_gaussians),s_tmp,t_tmp);
            if P_tmp < Poi
                qmin = q;
                q = (qmax + qmin)/2;                
            elseif P_tmp > Poi
                qmax = q;
                q = (qmax + qmin)/2;
            end
        end
        t(goi) = q;
        min_thresholds(1,:) = s;
        max_thresholds(1,:) = t;
        
        % determine thresholds for intermediate components
        for aix = 2:num_components-1
        coi = aix;
        goi = aix;
        s = -inf * ones(1, num_gaussians);
        for i = 1:coi-1
            s(i) = max_thresholds(i,i);
        end
        t =  inf * ones(1, num_gaussians);
        q = 0;
        qmin = -6;
        qmax = 6;
        Poi = P(coi);
        s_tmp = s;
        t_tmp = t; t_tmp(goi) = q;
        P_tmp = qsimvn(20,eye(num_gaussians),s_tmp,t_tmp);
        while abs(Poi - P_tmp) > P_tol
            t_tmp(goi) = q;
            P_tmp = qsimvn(20,eye(num_gaussians),s_tmp,t_tmp);
            if P_tmp < Poi
                qmin = q;
                q = (qmax + qmin)/2;                
            elseif P_tmp > Poi
                qmax = q;
                q = (qmax + qmin)/2;
            end
        end
        t(goi) = q;
        min_thresholds(coi,:) = s;
        max_thresholds(coi,:) = t;
        end
                
        % determine thresholds for last component
        coi = num_components;
        goi = num_components;
        s = -inf * ones(1, num_gaussians);
        for i = 1:coi-1
            s(i) = max_thresholds(i,i);
        end
        t =  inf * ones(1, num_gaussians);
        min_thresholds(coi,:) = s;
        max_thresholds(coi,:) = t;
end
