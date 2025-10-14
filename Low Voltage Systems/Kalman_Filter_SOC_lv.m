classdef Kalman_Filter_SOC_lv < matlab.System
    
    % Public, tunable properties
    properties
       qgain_1=1; 
       qgain_2=1;
       qgain_soc=1;
       N = 10;
    end

    % Pre-computed constants or internal states
    properties (Access = private)
        % Battery Lookup Table
        bat_lookup = load('battery_lookup_lv.mat').total; % Column 1 is SOC, 2 is OCV, 3 is dVdS (gradient)
        
        %State Space Matrices
        a
        b
        c

        %Process Covariance Update
        d
        d_window
        D
        
        %Output Equation
        g

        %Kalman Gain
        k
        
        %Battery Values
        r0
        r_ct
        r_dif
        c_ct
        c_dif
        
        %Battery Table Functions
        ke
        ocv
        
        q_nom = 3*3.5; %Capacity of three parallel cells, 3500mAh

        % Error covariance matrices
        q_fixed   % Process noise
        r %sensor noise
        
        % Initial conditions
        x
        p
        
        %Kalman Equations
        predX
        predP
        updateX
        updateK
        updateP

        %Time and Counter
        t
        j
       
    end

    methods (Access = protected)
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            obj.x = [0;0;0]; % Initial state, SoC,Vct,Vdif
            obj.p = [0 0 0;0 0 0;0 0 0]; % Initial error covariance
            obj.k = [1 1 1];
            obj.j=1;

            %Process Covariance Update
            obj.d_window = zeros(1,obj.N);

            obj.g = @(x, I) obj.ocv(x(1)) - x(2) - x(3) - obj.r0*I;
            
            %Update Equations
            obj.updateK = @(p,c,r) (c*p*c' + r)\p*c';      % K[k]   = PC^T/(CPC^T + R)
            obj.updateX = @(x,k,z,g) x + k*(z-g);        % x[n|n] = x[n|n-1] + K[n](z[n] - Cx[n|n-1]-Du[n]), old: obj.updateX = @(x,k,z,g) x + k*(z-c*x-Du);
            obj.updateP = @(k,c,p) (eye(size(p)) - k*c)*p; % P[n|n] = (I-K[n]C)P[n|n-1]
            
            %Sensor Error Matrices Initialization
            obj.r = 1e-4; %some value based on measurement or something, currently from paper
            obj.q_fixed = [obj.qgain_soc*1000*obj.r 0 0;0 obj.qgain_1*0.1*obj.r 0;...
                0 0 obj.qgain_2*0.01*obj.r];

            %Second Order Battery Model Numbers
            obj.r0 = 0.038; %lookup tables for all these constants based on temperature maybe they exist already
            obj.r_ct = 0;
            obj.r_dif = 0; 
            obj.c_ct = 0;
            obj.c_dif = 0;

            %Battery Table SoC, OCV Functions
            obj.ke = @(V) interp1(obj.bat_lookup(:,2), obj.bat_lookup(:,3), V, 'pchip','extrap'); %place/slope? on curve of OCV/SOC interpolate the slope
            obj.ocv = @(soc) interp1(obj.bat_lookup(:,1), obj.bat_lookup(:,2), soc, 'pchip','extrap');
            
            % Time/Prediction updates;
            obj.predX = @(a,x,b,u) a*x*a' + b*u;    % maybe not a'?? x[n+1|n] = Ax[n] + Bu[n]
            obj.predP = @(a,p,q) a*p*a' + q;        % P[n|n]   = AP[n|n]A^T + BQB^T
            obj.t = 0;
        end

        function [v1,v2,SoC, e_ekf, v, V_ocv] = stepImpl(obj, V, i,T,t_new)
            % Update Matrices
            obj.a = [0 0 1;exp(-(t_new-obj.t)/(obj.r_ct*obj.c_ct)) 0 0;...
                0 exp(-(t_new-obj.t)/(obj.r_dif*obj.c_dif))];
            obj.b = [-(t_new-obj.t)/obj.q_nom; obj.r_ct*(1-exp(-(t_new-obj.t)/(obj.r_ct*obj.c_ct)));...
                obj.r_dif*(1-exp(-(t_new-obj.t)/(obj.r_dif*obj.c_dif)))];
            obj.c = [obj.ke(V + obj.r0*i + sum([0 1 1]*obj.x,'All')) -1 -1];
            

            obj.t = t_new;
            
            %Predict
            obj.x = obj.predX(obj.a,obj.x,obj.b,i);
            obj.p = obj.predP(obj.a,obj.p,q);
            
            %Update
            obj.k = obj.updateK(obj.p, obj.c, obj.r);
            obj.x = obj.updateX(obj.x, obj.k, V, obj.g(obj.x,i)); %obj.x = obj.updateX(obj.x, obj.k, V, obj.c, obj.vocv0 + obj.d*i);
            obj.p = obj.updateP(obj.k, obj.c, obj.p);
          
            
            if obj.c*obj.x < 0
                obj.x(1) = 0;
            end

            %General Output
            obj.j=obj.j+1;
            v1= obj.x(1);
            %cov_v = sqrt(obj.c_v*obj.p*obj.c_v');
            v2 = obj.x(2);
            SoC = obj.(3);
            V_ocv = obj.ocv(SoC);
            v = obj.g(obj.x, i);
            obj.d = V - v;
            e_ekf = obj.d;

            %Process Covariance Update
            obj.d_window(1) = [];
            obj.d_window = [obj.d_window obj.d];
            obj.D = (1/obj.N)*sum(obj.d_window.*obj.d_window);
            obj.q_fixed = obj.k*obj.D*obj.k';

            
          
        end

        function resetImpl(~)
            % Initialize / reset internal properties
            %obj.x = [0;0]; % Initial state, v,a
            %obj.p = [0 0;0 0]; % Initial error covariance
            %obj.i=1;
        end
    end
end
