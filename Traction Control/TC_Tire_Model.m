global k0;
%global a0;

global C_x;
global B_x;
global E_x;
global C_y;
global B_y;
global E_y;
global D_x;
global D_y;
global fz_x_points;
global fz_y_points;
global fx_d_list;
global fy_d_list; 
global eps_k;
eps_k = 1e-12;

C_x = fx_fixed_params(1);
B_x = fx_fixed_params(2);
E_x = fx_fixed_params(3);

C_y = fy_fixed_params(1);
B_y = fy_fixed_params(2);
E_y = fy_fixed_params(3);

%READ, MUST LOAD VARIABLES FROM MF_var_coeff and MF_fx_coeff, Benjamin
%version to work
%READ: For testing with previous tire data
%tire_data = load("pick your file"); % 18in R20
SA = tire_data.SA;
FZ = tire_data.FZ;
FY = tire_data.FY;
IA = tire_data.IA;
FX = tire_data.FX;
SR = tire_data.SL;
P = tire_data.P;

Fy_list = zeros(1,numel(SA));
Fx_list = zeros(1,numel(SA));

Fsy_list = zeros(1,numel(SA)); 
Fay_list = zeros(1,numel(SA)); 
f0y_list = zeros(1,numel(SA)); 
for n = 1:numel(FZ)
    fz = FZ(n);
        
    [Fx,Fy,F_sy,F_ay,f0y,f0y2] = compute_forces(SR(n),deg2rad(SA(n)),fz);
    Fy_list(n) = Fy;
    Fx_list(n) = Fx;
    Fsy_list(n) = F_sy; 
    Fay_list(n) = F_ay;   
    f0y_list(n) = f0y;  

end
FY_lookup = [SR SA Fy_list'];
FX_lookup = [SR SA Fx_list'];

%{
READ: For varying z, a, and k on their own for the model
fz_values = (0:500:5000);
fz_names = {};
for m = 1:numel(fz_values)
    fz_names{m} = fz_values(m);
end
FY_lookup = struct();  
FX_lookup = struct();% start empty
for k = 1:numel(fz_names)
    FY_lookup.fz_names{k} = [];   % initialize with empty
    FX_lookup.fz_names{k} = [];
end
count = 0;
for n = 1:numel(fz_values)
    fz = fz_values(n);
    disp('here');
    disp(count);
    
    k_vary = (-1:0.01:1);
    a_vary = (-30:0.2:30);
    [a_list, k_list] = meshgrid(a_vary, k_vary);
    a_list = a_list(:);
    k_list = k_list(:);
    numel(a_list);
    Fy_list = zeros(1,numel(a_list));
    Fx_list = zeros(1,numel(a_list));
   
    Fsy_list = zeros(1,numel(a_list)); % ADD THIS
    Fay_list = zeros(1,numel(a_list)); % ADD THIS
    f0y_list = zeros(1,numel(a_list)); % ADD THIS
    for i=1:numel(a_list)
        
        [Fx,Fy,F_sy,F_ay,f0y,f0y2] = compute_forces(k_list(i),deg2rad(a_list(i)),fz);
        Fy_list(i) = Fy;
        Fx_list(i) = Fx;
        Fsy_list(i) = F_sy;   % ADD THIS
        Fay_list(i) = F_ay;   % ADD THIS
        f0y_list(i) = f0y;   % ADD THIS
        
    end
    disp('here2')
    FY_lookup.fz_names{n} = [k_list a_list Fy_list'];
    FX_lookup.fz_names{n} = [k_list a_list Fx_list'];
    count=count+1;
end
%}
function [f_x,f_y,f_sy,f_ay,f0y,f0y2] = compute_forces(k,a,fz)
    global D_x;
    global D_y;

    [D_x,D_y] = compute_ds(fz);
    [f_ax,f_ay,f0y] = run_adhere_forces(k,a);
    [f_sx,f_sy,f0y2] = run_slip_forces(k,a);
    
    f_x = f_sx + f_ax;
    f_y = f_sy + f_ay;
end

function [f_ax,f_ay,f0y] = run_adhere_forces(k,a)
    [sig_xo,sig_yo] = compute_sigmaso();
    
    expression1 = ((k/sig_xo)*cos(a))^2 + (sin(a)/sig_yo)^2;
    expression2 = (1-k)^2*(cos(a))^2;
    if expression1 < expression2
        [sig_x,sig_y] = compute_sigmas(k,a);
        psi = compute_psi(sig_x,sig_y,sig_xo,sig_yo);
        [k0,a0] = compute_pure_slips(k,a);
        [f0x,f0y] = compute_pacejka_forces(k0,a0,1);
        [sig_x0,sig_y0] = compute_sigmas0(k0,a0);
        
        psi_x0 = compute_psi(sig_x0,0,sig_xo,sig_yo);
        psi_y0 = compute_psi(0,sig_y0,sig_xo,sig_yo);
        [f_ax,f_ay] = compute_adhere_forces(f0x,f0y,psi,psi_x0,psi_y0);
    else
        f_ax = 0;
        f_ay = 0;
        f0y = 0;
    end
end

function [f_sx,f_sy,f0y] = run_slip_forces(k,a)
    
    [k0_vel,a0_vel] = compute_slipvel_pure_slips(k,a);
    [f0x,f0y] = compute_pacejka_forces(k0_vel,a0_vel,0);
    [sig_xo,sig_yo] = compute_sigmaso();
    B = compute_B(k,a);
    sig_x = 0;
    sig_y = 0;
    Bo = 0;
    psi = 1;
    expression1 = ((k/sig_xo)*cos(a))^2 + (sin(a)/sig_yo)^2;
    expression2 = (1-k)^2*(cos(a))^2;
    if expression1 < expression2
        
        [sig_x,sig_y] = compute_sigmas(k,a);
        psi = compute_psi(sig_x,sig_y,sig_xo,sig_yo);
        
        psi_x = compute_psi(sig_x,0,sig_xo,sig_yo);
        psi_y = compute_psi(0,sig_y,sig_xo,sig_yo);
        Bo = compute_Bo(psi_x,psi_y);
    end
    expression1 = sqrt((k*cos(a))^2+(sin(a))^2);
    denom_check = 1-expression1*sign(k);
    if denom_check > 0 && expression1 < sig_xo*abs(denom_check) 
        sig_x0_vel = compute_sig_x0_vel(k,a);
        psi_x0_vel = compute_psi(sig_x0_vel,0,sig_xo,sig_yo);
        
        theta_x = compute_theta_x(psi_x0_vel);
    else
        theta_x = 1;
        psi_x0_vel = 1;
    end
    
    denom_check2 = 1-expression1^2;
    if denom_check2 > 0 && expression1 < sig_yo*sqrt(denom_check2)
        sig_y0_vel = compute_sig_y0_vel(k,a);
        psi_y0_vel = compute_psi(0,sig_y0_vel,sig_xo,sig_yo);

        theta_y = compute_theta_y(psi_y0_vel);
    else
        theta_y = 1;
        psi_y0_vel = 1;
    end
    
    [Gamma_x,Gamma_y] = compute_gammas(sig_x,sig_y,sig_xo,sig_yo,B,Bo,psi,psi_x0_vel,psi_y0_vel);

    f0_sx = theta_x*f0x*Gamma_x;
    f0_sy = theta_y*f0y*Gamma_y;
    Bf = compute_Bf(k,a,f0_sx,f0_sy);
   
    [f_sx,f_sy] = compute_slip_forces(f0_sx,f0_sy,Bf);

end

function [D_x, D_y] = compute_ds(fz)
    global fz_x_points;
    global fz_y_points;
    global fx_d_list;
    global fy_d_list;
    
    D_x = interp1(fz_x_points, fx_d_list, fz,'makima','extrap');
    D_y = interp1(fz_y_points, fy_d_list, fz, 'makima','extrap');
    
    if isnan(D_x) || isinf(D_x)
        D_x
    elseif isnan(D_y) || isinf(D_y)
        D_y
    end
end


function psi = compute_psi(sig_x,sig_y,sig_x0,sig_y0)
    global B_x;
    global C_x;
    global B_y;
    global C_y;
    global k0;
    global D_x;
    global D_y;
    global eps_k;
 
    psi = sqrt((sig_x/sig_x0)^2+(sig_y/sig_y0)^2);
    
end

function [sig_x, sig_y] = compute_sigmas(k,a) 
    global eps_k;
    if abs(1-k) < eps_k
        den = eps_k;
    else
        den = 1-k;
    end
    sig_x = k/den;
    sig_y = tan(a)/den;
end

function [sig_x, sig_y] = compute_sigmaso()
    global B_x;
    global C_x;
    global B_y;
    global C_y;
    global D_x;
    global D_y;
    
    sig_x = 3/(B_x*C_x);
    sig_y = D_y*((2/(B_x*C_x*D_x))+(1*pi/(B_y*C_y*D_y*180)));
end

function [k0,a0] = compute_pure_slips(k,a)
    k0 = k;
    a0 = atan2(sin(a),(1-k)*cos(a)); 
    
end

function [f0x,f0y] = compute_pacejka_forces(k,a,check)

    global B_x;
    global C_x;
    global E_x;
    global B_y;
    global C_y;
    global E_y;
    global D_x;
    global D_y;
    a = rad2deg(a);
    
    if check == 1
        if k < 0
            k_a = (-k)/(1-2*k);
        else
            k_a = k;
        end  
    else
        if k<0
            k_a = -k;
        else
            k_a = k;
        end
    end
    
    if a < 0
        a_a = -a;
    else
        a_a = a;
    end
    
    f0x = D_x*sin(C_x*atan(B_x*(1-E_x)*k_a+E_x*atan(B_x*k_a)))*sign(-k);
    f0y = D_y*sin(C_y*atan(B_y*(1-E_y)*a_a+E_y*atan(B_y*a_a)))*sign(a);

end

function [sig_x0, sig_y0] = compute_sigmas0(k,a)
    global eps_k;
    if abs(1-k) < eps_k
        den = eps_k;
    else
        den = 1-k;
    end
    sig_x0 = k/den;
    sig_y0 = tan(a)/den;

end

function [f_ax, f_ay] = compute_adhere_forces(f0x,f0y,psi,psi_x0,psi_y0)
    global eps_k;
    denom_x = 3*(1-psi_x0)^2+psi_x0*(3-2*psi_x0);
    denom_y = 3*(1-psi_y0)^2+psi_y0*(3-2*psi_y0);
    if abs(denom_x) < eps_k
        f_ax = 0;
    else
        f_ax = ((3*(1-psi)^2)*f0x)/denom_x;
    end
    
    if abs(denom_y) < eps_k
        f_ay = 0;
    else
        f_ay = ((3*(1-psi)^2)*f0y)/denom_y;
    end
    
end

function [k0_vel,a0_vel] = compute_slipvel_pure_slips(k,a)
    s = hypot(k*cos(a), sin(a));
    s = max(s, 1e-12); 
    k0_vel = s*sign(k); %v/v0 = 1 
    a0_vel = asin(sat(s*sign(a),-1,1)); %v/v0 = 1
end

function B = compute_B(k,a)
    B = atan2(tan(a),k);
end

function Bo = compute_Bo(psi_x,psi_y)
    Bo = atan2(psi_y,psi_x);
end

function sig_x0_vel = compute_sig_x0_vel(k,a)

    expression1 = sqrt((k*cos(a))^2+(sin(a))^2)*sign(k);
    expression2 = 1-sqrt((k*cos(a))^2+(sin(a))^2)*sign(k);
    sig_x0_vel = expression1/expression2;

end

function sig_y0_vel = compute_sig_y0_vel(k,a)

    expression1 = sqrt((k*cos(a))^2+(sin(a))^2)*sign(a);
    expression2 = sqrt(1-((k*cos(a))^2+(sin(a))^2));
    sig_y0_vel = expression1/expression2;

end

function theta_x = compute_theta_x(psi)
    expression1 = psi*(3-2*psi);
    expression2 = 3*(1-psi)^2+psi*(3-2*psi);
    theta_x = expression1/expression2;
end

function theta_y = compute_theta_y(psi) 
    expression1 = psi*(3-2*psi);
    expression2 = 3*(1-psi)^2+psi*(3-2*psi);
    theta_y = expression1/expression2;
end

function [Gamma_x,Gamma_y] = compute_gammas(sig_x,sig_y,sig_xo,sig_yo,B,Bo,psi,psi_x0_vel,psi_y0_vel)
    if psi < 1
        if psi_x0_vel < 1
            expression1 = ((3-2*psi)/(3-2*psi_x0_vel))...
                *(sig_xo*((1/sig_xo)*abs(cos(B))*cos(Bo)...
                +(1/sig_yo)*abs(sin(B))*sin(Bo)))^2;
            expression2 = (sqrt((1+sig_x)^2+sig_y^2)-sqrt(sig_x^2+sig_y^2)*sign(sig_x))^2;
            Gamma_x = expression1*expression2;
        else
            Gamma_x = psi^2*(3-2*psi);
        end

        if psi_y0_vel < 1
            expression1 = ((3-2*psi)/(3-2*psi_y0_vel))...
                *(sig_yo*((1/sig_xo)*abs(cos(B))*cos(Bo)...
                +(1/sig_yo)*abs(sin(B))*sin(Bo)))^2;
            expression2 = (((1+sig_x)^2+sig_y^2)-(sig_x^2+sig_y^2));
            Gamma_y = expression1*expression2;
        else
            Gamma_y = psi^2*(3-2*psi);
        end
    else
        if psi_x0_vel < 1
            Gamma_x = 1/(psi_x0_vel^2*(3-2*psi_x0_vel));
        else
            Gamma_x = 1;
        end

        if psi_y0_vel < 1
            Gamma_y = 1/(psi_y0_vel^2*(3-2*psi_y0_vel));
        else
            Gamma_y = 1;
        end
    end
end

function Bf = compute_Bf(k,a,f0_sx,f0_sy)

    global eps_k;

    vsx = k*cos(a);  vsy = sin(a);
    if abs(vsy) < eps_k
        Bf = 0;
    elseif abs(vsx) < eps_k
        Bf = (pi/2)*sign(a);
    else
        Bf = atan2(f0_sx*vsy, f0_sy*vsx);
    end
 
end


function [f_sx,f_sy] = compute_slip_forces(f0_sx,f0_sy,Bf)
    f_sx = f0_sx*abs(cos(Bf));
    f_sy = f0_sy*abs(sin(Bf));
end

function y = sat(u,minV,maxV)
    y = min(maxV,max(minV,u));
end

