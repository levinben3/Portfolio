classdef Kalman_Filter_SOC < matlab.System
    % untitled3 Add summary here
    %
    % This template includes the minimum set of functions required
    % to define a System object.

    % Public, tunable properties
    properties
       qgain_1=1; 
       qgain_2=1;
       qgain_soc=1;
    end

    % Pre-computed constants or internal states
    properties (Access = private)
        % State-space matrices
        bat_lookup = load('battery_lookup.mat'); % Column 1 is SOC, 2 is OCV, 3 is dVdS (gradient)
        a; 
        b; 
        c;
        d;
        k=1;
        r0
        r1
        r2
        c1
        c2
        vocv0
        ke
        S_c = 3*3.5; %3*capacity of three parallel cells, 3500mAh
        % Error covariance matrices
        q_fixed   % Process noise
        r %sensor noise
        
        % Initial conditions
        x
        p
        y
        predX
        predP
        updateX
        updateK
        updateP
        t
       
    end

    methods (Access = protected)
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            obj.x = [0;0;0]; % Initial state, V1,V2,SoC
            obj.p = [0 0 0;0 0 0;0 0 0]; % Initial error covariance
            obj.y = 0;
            obj.j=1;
            obj.r = 8.432e-4; %some value from a paper, seems it would be small
            
            obj.updateK = @(p,c,r) (c*p*c' + r)\p*c';      % K[k]   = PC^T/(CPC^T + R)
            obj.updateX = @(x,k,z,c,Du) x + k*(z-c*x-Du);        % x[n|n] = x[n|n-1] + K[n](z[n] - Cx[n|n-1]-Du[n])
            obj.updateP = @(k,c,p) (eye(size(p)) - k*c)*p; % P[n|n] = (I-K[n]C)P[n|n-1]
            obj.r = 0; %some value based on measurement or something
            obj.r0 = @(T) 0; %lookup tables for all these constants based on temperature maybe they exist already
            obj.r1 = @(T) 0;
            obj.r2 = @(T) 0; 
            obj.c1 = @(T) 0;
            obj.c2 = @(T) 0;
            obj.ke = @(V) interp1(obj.bat_lookup(:,1), obj.bat_lookup(:,2), V, 'pchip','extrap'); %place/slope? on curve of OCV/SOC interpolate the slope
            % Time updates
            obj.predX = @(a,x,b,u) a*x + b*u;            % x[n+1|n] = Ax[n] + Bu[n]
            obj.predP = @(a,p,q) a*p*a' + q;        % P[n|n]   = AP[n|n]A^T + BQB^T
            obj.t = 0;
        end

        function [v1,v2,SoC] = stepImpl(obj, V, i,T,t_new)
            % Implement algorithm. Calculate y as a function of input u and
            % internal states.
            obj.a = [exp(-(t_new-obj.t)/(obj.r1(T)*obj.c1(T))) 0 0;...
            0 exp(-(t_new-obj.t)/(obj.r2(T)*obj.c2(T))) 0;0 0 1];
            obj.b = [obj.r1(T)*(1-exp(-(t_new-obj.t)/(obj.r1(T)*obj.c1(T))));...
            obj.r2(T)*(1-exp(-(t_new-obj.t)/(obj.r2(T)*obj.c2(T)))); -(t_new-obj.t)/obj.S_c];
            obj.d = -obj.r0(T);
            obj.c = [-1 -1 obj.ke(V - obj.d*i + [1 1 0]*obj.x)];
            obj.q_fixed = [obj.qgain_1*0.1*obj.r 0 0;0 obj.qgain_2*0.01*obj.r 0;...
            0 0 obj.qgain_soc*1000*obj.r];
            q=obj.qgain*obj.q_fixed;

            obj.t = t_new;
            
            %Predict
            obj.x = obj.predX(obj.a,obj.x,obj.b,i);
            obj.p = obj.predP(obj.a,obj.p,q);
           
            %Update
            obj.k = obj.updateK(obj.p, obj.c, obj.r);
            obj.x = obj.updateX(obj.x, obj.k, V, obj.c, obj.vocv0 + obj.d*i);
            obj.p = obj.updateP(obj.k, obj.c, obj.p);
          
        
            if obj.c*obj.x < 0
                obj.x(1) = 0;
            end

            obj.j=obj.j+1;
            v1= obj.x(1);
            %cov_v = sqrt(obj.c_v*obj.p*obj.c_v');
            v2 = obj.x(2);
            SoC = obj.(3);
          
        end

        function resetImpl(~)
            % Initialize / reset internal properties
            %obj.x = [0;0]; % Initial state, v,a
            %obj.p = [0 0;0 0]; % Initial error covariance
            %obj.i=1;
        end
    end
end
