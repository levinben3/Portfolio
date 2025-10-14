classdef Kalman_Filter_SOC_hv < matlab.System
    
    % Public, tunable properties
    properties
       qgain_1=1; 
       qgain_2=1;
       qgain_soc=1;
       qgain=1;
    end

    % Pre-computed constants or internal states
    properties (Access = private)
        % Battery Lookup Table
        bat_lookup = load('battery_lookup_slopes.mat').total; % Column 1 is SOC, 2 is OCV, 3 is dVdS (gradient)
        
        %State Space Matrices
        a
        b
        c
        
        %Process Covariance Update
        d
        d_window
        D

        %Output equation
        g
        
        %Parameter Matrices
        c_theta

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

        q_nom = 8*2.6; %Capacity of eight parallel cells, 2600mAh
        
        % Error covariance matrices
        q_fixed   % Process noise
        q_theta_fixed
        r %sensor noise
        r_theta
        
        % Initial conditions
        dx_dtheta
        I_prev
        x
        theta
        p
        p_theta
        
        
        %Kalman Equations
        predX
        predP
        predTheta
        predPTheta
        updateX
        updateK
        updateP

        %Time and Counter
        t
        j
       
    end

    methods (Access = protected)
        function setupImpl(obj)
            
            
            % Initialize Vectors
            obj.x = [0;0;0]; % Initial state, SoC,Vct,Vdif
            obj.theta = [obj.r0;obj.q_nom]; % Parameters Vector R_0 Q_nom
            obj.p = [0 0 0;0 0 0;0 0 0]; % Initial error covariance
            obj.p_theta = [0 0;0 0];
            obj.k = [1 1 1];
            
            obj.j=1;

            %Process Covariance Update
            obj.d_window = zeros(1,obj.N);
            
            obj.I_prev = 0;
            obj.dx_dtheta = zeros(3,2);
            
            obj.g = @(x, I, theta) obj.ocv(x(1)) - x(2) - x(3) - theta(1)*I;

            %Update Equations
            obj.updateK = @(p,c,r) (c*p*c' + r)\p*c';      % K[k]   = PC^T/(CPC^T + R)
            %obj.updateX = @(x,k,z,c,Du) x + k*(z-c*x-Du);        % x[n|n] = x[n|n-1] + K[n](z[n] - Cx[n|n-1]-Du[n])
            obj.updateX = @(x,k,z,g) x + k*(z-g);        % MAYBE x[n|n] = x[n|n-1] + K[n](z[n] - Cx[n|n-1]-Du[n])
            obj.updateP = @(k,c,p) (eye(size(p)) - k*c)*p; % P[n|n] = (I-K[n]C)P[n|n-1]
            
            %Sensor Error Matrices Initialization
            obj.r = 1e-4; %from paper, think about, probably small
            obj.r_theta = 10; %from paper
            obj.q_fixed = [obj.qgain_soc*1000*obj.r 0 0;0 obj.qgain_1*0.1*obj.r 0;...
                0 0 obj.qgain_2*0.01*obj.r];
            
            %Second Order Battery Model Numbers
            obj.r0 = 0.04; %Find this value
            obj.r_ct = @(T) 0;%lookup tables for all these constants based on temperature maybe they exist already
            obj.r_dif = @(T) 0; 
            obj.c_ct = @(T) 0;
            obj.c_dif = @(T) 0;

            %Battery Table SoC, OCV Functions
            obj.ke = @(V) interp1(obj.bat_lookup(:,1), obj.bat_lookup(:,2), V, 'pchip','extrap'); %place/slope? on curve of OCV/SOC interpolate the slope
            obj.ocv = @(soc) interp1(obj.bat_lookup(:,1), obj.bat_lookup(:,2), soc, 'pchip','extrap');
            
            % Time/Prediction updates;
            obj.predX = @(a,x,b,u) a*x*a' + b*u;    % maybe not a'?? x[n+1|n] = Ax[n] + Bu[n]
            obj.predP = @(a,p,q) a*p*a' + q;        % P[n|n]   = AP[n|n]A^T + BQB^T
            obj.predTheta = @(theta) theta;
            obj.predPTheta = @(p,q) p + q;
            obj.t = 0;
        end

        function [v1,v2,SoC, e_ekf, v, V_ocv] = stepImpl(obj, V, i,T,t_new)
            % Implement algorithm. Calculate y as a function of input u and
            % internal states.
            obj.a = [0 0 1;exp(-(t_new-obj.t)/(obj.r_ct(T)*obj.c_ct(T))) 0 0;...
                0 exp(-(t_new-obj.t)/(obj.r_dif(T)*obj.c_dif(T)))];
            obj.b = [-(t_new-obj.t)/obj.theta(2); obj.r_ct(T)*(1-exp(-(t_new-obj.t)/(obj.r_ct(T)*obj.c_ct(T))));...
                obj.r_dif(T)*(1-exp(-(t_new-obj.t)/(obj.r_dif(T)*obj.c_dif(T))))];
            obj.c = [obj.ke(V + obj.theta(1)*i + sum([0 1 1]*obj.x,'All')) -1 -1];
            obj.c_theta = [-i 0] + obj.c*([0 ((t_new-obj.t)*obj.I_prev)/(obj.theta(2))^2;0 0;0 0]+obj.a*obj.dx_dtheta);
            obj.dx_dtheta = [0 ((t_new-obj.t)*obj.I_prev)/(obj.theta(2))^2;0 0;0 0] + obj.a*obj.dx_dtheta;
            
            obj.I_prev = i;
            obj.q_theta_fixed = obj.qgain*[1e-6 0;0 1e-4];
            obj.t = t_new;
            
            %Predict
            obj.x = obj.predX(obj.a,obj.x,obj.b,i);
            obj.p = obj.predP(obj.a,obj.p,q);
            
            obj.theta = obj.predTheta(obj.theta);
            obj.p_theta = obj.predPTheta(obj.p_theta, obj.qgain*obj.q_theta_fixed);
           
            %Update
            obj.k = obj.updateK(obj.p, obj.c, obj.r);
            obj.x = obj.updateX(obj.x, obj.k, V, obj.g(obj.x, I, obj.theta));
            obj.p = obj.updateP(obj.k, obj.c, obj.p);

            obj.k = obj.updateK(obj.p_theta, obj.c_theta, obj.r_theta);
            obj.theta = obj.updateX(obj.theta, obj.k, V, obj.g(obj.x, I, obj.theta));
            obj.p_theta = obj.updateP(obj.k, obj.c_theta, obj.p_theta);
          
        
            if obj.c*obj.x < 0
                obj.x(1) = 0;
            end

            if obj.c_theta*obj.theta < 0
                obj.theta(1) = 0;
            end

            obj.j=obj.j+1;
            v1= obj.x(1);
            %cov_v = sqrt(obj.c_v*obj.p*obj.c_v');
            v2 = obj.x(2);
            SoC = obj.(3);
            V_ocv = obj.ocv(SoC);
            v = obj.g(obj.x, i, obj.theta);
            e_ekf = V - v;

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
