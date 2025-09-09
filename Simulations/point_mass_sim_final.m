%% Regulations
fsae_ev.power_limit = 80e3;
fsae_hybrid.power_limit = inf;

%% Tracks
icahn_loop.description = "Icahn Stadium Test Loop";
icahn_loop.location = "New York, NY";
icahn_loop.air_density = 1.293;
icahn_loop.segments.lengths   = [60 11.75*pi 20 15.25*pi/2 10 18.75*pi/2 40 11.75*pi/2 19 21.75*pi/2 14 3.25*pi/2 5 6.75*pi/2];
icahn_loop.segments.curvature = [0  -1/11.75 0  1/15.25    0  -1/18.75   0  -1/11.75   0  -1/21.75   0  1/3.25    0 -1/6.75];
icahn_loop.segments.limits    = Inf * ones(1, length(icahn_loop.segments.lengths));

accel_event.description = "FSAE Acceleration Event";
accel_event.location = "Brooklyn, MI";
accel_event.air_density = 1.293;
accel_event.segments.lengths   = [10 0.3 75  75  10];
accel_event.segments.curvature = [0  0   0   0   0];
accel_event.segments.limits    = [0.1  Inf Inf Inf 0.1];

skidpad_event.description = "FSAE Skidpad Event";
skidpad_event.location = "Brookyln, MI";
skidpad_event.air_density = 1.293;
R = 18.25/2; % min 15.25/2 + track/2, max 21.25/2 - track/2
skidpad_event.segments.lengths   = [R R   2*pi*R  2*pi*R 2*pi*R 2*pi*R R];
skidpad_event.segments.curvature = [0 0   1/R     1/R    -1/R   -1/R   0];
skidpad_event.segments.limits    = [0.1 Inf Inf     Inf    Inf    Inf  Inf];
% skidpad_event.segments.lengths   = [2*pi*R];
% skidpad_event.segments.curvature = [1/R];
% skidpad_event.segments.limits    = [Inf];

endur_event.description = "FSAE Endurance Event Michigan 2024";
endur_event.location = "Brookln, MI";
endur_event.air_density = 1.293;
endur_event.segments.lengths = [7 5*pi/6 2*pi/2.8 3 2*pi/3 4 1.5*pi/6 1.5 1*pi/6 1.4*pi/4 1.2*pi/3 0.7*pi/5 1.5 1.5*pi/3 1*pi/2.5 5*pi/10 5*pi/18 7 ...
1.5*pi/10 1*pi/2.2 1 1.5*pi/2.8 3 1*pi/1.2 0.9 1.5*pi/2.5 0.7 0.5*pi/1.8 5*pi/12 2 6.5*pi/3.5 2*pi/3.5 1.5 1*pi/3.5 3*pi/6.5 2 ...
0.7*pi/5 0.5*pi/4.6 5 3*pi/9 3*pi/5 1*pi/3.5 1*pi/2.5 1*pi/2.5 1*pi/2.5 1.5*pi/3.25 2.5 1*pi/2.3 1.7*pi/1.3 1.3*pi/2.6 3.5 0.5*pi/4 ...
1*pi/2.3 1*pi/2.5 1*pi/2.3 1 0.5*pi/2.4 1.3 0.5*pi/2.3 1.5 0.8*pi/2.1 1 0.8*pi/2.1 0.8 0.9*pi/1.4 2.1 1.5*pi/3 1.5 1*pi/2.7 1.7 ...
1.2*pi/4.5 2 1*pi/3.4 1.5454 2.4123];
endur_event.segments.curvature = [0 1/5 -1/2 0 1/2 0 -1/1.5 0 1/1 -1/1.2 1/1.2 -1/0.7 0 -1/1.5 1/1 -1/5 1/5 0 ...
-1/1.5 1/1 0 -1/1.5 0 -1/1 0 1/1.5 0 -1/0.5 1/5 0 -1/7 1/2 0 -1/1 1/3 0 ...
-1/0.7 1/0.5 0 1/3 -1/3 1/1 -1/1 1/1 -1/1 1/2 0 -1/1 1/1.7 -1/1.3 0 -1/0.5 ...
1/1 -1/1 1/1 0 -1/0.5 0 1/0.5 0 -1/0.8 0 1/0.8 0 -1/0.9 0 -1/1.5 0 -1/1 0 ...
1.2/1 0 -1/1 0 1/3];
factor = 1000 / sum(endur_event.segments.lengths);
endur_event.segments.lengths = factor * endur_event.segments.lengths;
endur_event.segments.curvature = 1/factor * endur_event.segments.curvature;
endur_event.segments.limits = Inf * ones(1, length(endur_event.segments.lengths));

%% Motors
emrax208.max_rpm = 7000;
emrax228.max_rpm = 5000;
emrax268.max_rpm = 4500;
dhx_k40.max_rpm = 6000;

% speed-torque characteristics:
emrax208.LUT_WM_W = 0.10472*[0;4000;4500;5000;5500;6000;7000];
emrax208.LUT_WM_M = [140;135;132;128;122;115;0];

emrax228.LUT_WM_W = 0.10472*[0;2000;3000;4000;5000;5001];
emrax228.LUT_WM_M = [240;240;234;228;216;0];

emrax268.LUT_WM_W = 0.10472*[0;2000;4500;4501];
emrax268.LUT_WM_M = [500;500;450;0];

dhx_k40.LUT_WM_W = 0.10472*[0;6000;6001];
dhx_k40.LUT_WM_M = [80;80;0];

% efficiency maps
emrax208.LUT_WME_W = 0.10472*[500;1000;1500;2000;2500;3000;3500;4000;4500;5000];
emrax208.LUT_WME_M = [20;40;60;80;100;120;140];
emrax208.LUT_WME_E = [86 88 89 88 86 84 82;88 94 94 94 92 89 85;89 94.5 95.5 96 94.5 91 85;89 95.5 96 96 95 92 85;89 95.5 96 96 95 92 85;89 95 96 95.75 94.75 92 85;88 94.5 96 95.5 94.5 92 85;86 94 95 95.25 94.5 92 85;80 87 94 94.5 94.5 91 85;80 80 86 94 94 90 86]/100;

emrax228.LUT_WME_W = 0.10472*[500;1000;1500;2000;2500;3000;3500;4000;4500;5000]*5500/6000;
emrax228.LUT_WME_M = [20;40;60;80;100;120;140]*240/140;
emrax228.LUT_WME_E = [86 88 89 88 86 84 82;88 94 94 94 92 89 85;89 94.5 95.5 96 94.5 91 85;89 95.5 96 96 95 92 85;89 95.5 96 96 95 92 85;89 95 96 95.75 94.75 92 85;88 94.5 96 95.5 94.5 92 85;86 94 95 95.25 94.5 92 85;80 87 94 94.5 94.5 91 85;80 80 86 94 94 90 86]/100;

emrax268.LUT_WME_W = 0.10472*[500;1000;1500;2000;2500;3000;3500;4000;4500;5000]*4500/6000;
emrax268.LUT_WME_M = [20;40;60;80;100;120;140]*500/140;
emrax268.LUT_WME_E = [86 88 89 88 86 84 82;88 94 94 94 92 89 85;89 94.5 95.5 96 94.5 91 85;89 95.5 96 96 95 92 85;89 95.5 96 96 95 92 85;89 95 96 95.75 94.75 92 85;88 94.5 96 95.5 94.5 92 85;86 94 95 95.25 94.5 92 85;80 87 94 94.5 94.5 91 85;80 80 86 94 94 90 86]/100;

dhx_k40.LUT_WME_W = [0; 7000];
dhx_k40.LUT_WME_M = [0; 80];
dhx_k40.LUT_WME_E = [1 1; 1 1] * .93; % todo: fill in efficiency map

% hard power limit per motor
emrax208.P = 41.5e3;
emrax228.P = 80e3;
emrax268.P = 60e3;
dhx_k40.P = 40e3;

%% Cars
ev24.mass = 264 + 68;
ev24.cg = [738 0 255.41];
ev24.aero.cda = 0.428;
ev24.aero.cla = 0;
ev24.drive.motor = emrax208;
ev24.drive.ratio = 3.7;
ev24.drive.efficiency = 0.96;
ev24.drive.count = 1;
ev24.tires.mux = 1.5; % TODO: get mu from chris
ev24.tires.muy = 1.5;
ev24.tires.rolling_resistance = 0.015;
ev24.tires.radius = 0.2032;
ev24.hv.vmax = 302.4;
ev24.hv.vnom = 260;

ev25.mass = 242.5 + 68;
ev25.cg = [738 0 255.41];
ev25.aero.cda = 0.428;
ev25.aero.cla = 0;
ev25.drive.motor = emrax208;
ev25.drive.ratio = 4.9;
ev25.drive.efficiency = 0.96;
ev25.drive.count = 1;
ev25.tires = ev24.tires;
ev25.hv = ev24.hv;

% 2 x 268
ev26a.mass = 260 + 68;
ev26a.cg = ev25.cg;
ev26a.l = 1.530;
ev26a.aero.cda = 0.2288;
ev26a.aero.cla = 0.51;
ev26a.drive.motor = emrax268;
ev26a.drive.ratio = 1;
ev26a.drive.efficiency = 0.96;
ev26a.drive.count = 2;
ev26a.tires = ev25.tires;
ev26a.hv.vmax = 600;
ev26a.hv.vnom = 520;

% 1 x 208
ev26b.mass = 220 + 68;
ev26b.cg = ev25.cg;
ev26b.l = 1.530;
ev26b.aero.cda = 0.2288; % will
ev26b.aero.cla = 0.51;
ev26b.drive.motor = emrax208;
ev26b.drive.ratio = 4.9;
ev26b.drive.efficiency = 0.96;
ev26b.drive.count = 1;
ev26b.tires = ev25.tires; % chris
ev26b.hv.vmax = 300;
ev26b.hv.vnom = 260;

%% Accel Experiment

t = accel_event;
i=1;
sim26a = summarize_lapsim(fsae_ev, t, ev26a);
sim26b = summarize_lapsim(fsae_ev, t, ev26b);
max(sim26b.stats.pelectric)
max(sim26b.vv)
figure(9)
subplot(211);
plot(sim26b.vv);
subplot(212);
plot(sim26b.stats.ax);
accel = [sim26a sim26b];
figure(Name="Accel Event");
plot_lap(sim26a, sim26b,i);
accel_times = arrayfun(@(s) s.splits(3) - s.splits(2), accel)

%% Skidpad Experiment
t = skidpad_event;
i = i+2;
sim26a = summarize_lapsim(fsae_ev, t, ev26a);
sim26b = summarize_lapsim(fsae_ev, t, ev26b);
skidpad = [sim26a sim26b];
figure(Name="Skidpad Event");
plot_lap(sim26a, sim26b,i);
skidpad_times = arrayfun(@(s) s.splits(4) - s.splits(3), skidpad)
%skidpad_times = arrayfun(@(s) s.time, skidpad)

%% Icahn Experiment
t = icahn_loop;
i = i+2;
sim26a = summarize_lapsim(fsae_ev, t, ev26a);
sim26b = summarize_lapsim(fsae_ev, t, ev26b);
icahn = [sim26a sim26b];
figure(Name="Icahn Loop");
plot_lap(sim26a, sim26b, i);
icahn_times = arrayfun(@(s) s.time, icahn)

%% Endurance Experiment
t = endur_event;
i = i+2;
sim26a = summarize_lapsim(fsae_ev, t, ev26a);
sim26b = summarize_lapsim(fsae_ev, t, ev26b);
max(sim26b.stats.pelectric)
max(sim26b.vv)
endur = [sim26a sim26b];
figure(Name="Endurance");
plot_lap(sim26a, sim26b, i);
endur_times = arrayfun(@(s) s.time * 22, endur)
figure(10)
subplot(211)
plot(sim26b.vv)
subplot(212);
plot(sim26b.stats.ax);
figure(11)
plot(sim26b.stats.ax,sim26b.stats.ay,'*');
energya = sim26a.stats.avgelectric*endur_times(1)/3600;
energyb = sim26b.stats.avgelectric*endur_times(2)/3600;
%% Competition Points
score26a = lap2score(accel_times(1), skidpad_times(1), endur_times(1) * 0.8/22, endur_times(1), energya)
score26b = lap2score(accel_times(2), skidpad_times(2), endur_times(2) * 0.8/22, endur_times(2), energyb)

%% Plotting
function plot_lap(sim26a, sim26b, i)
    M = 4;
    N = 2;
    P = 1;
    figure(i)
    subplot(M,N,P); P = P + 1;
    plot(sim26a.xx, sim26a.vv * 3.6);
    xlabel("Distance (m)")
    ylabel("Velocity (km/h)")
    hold on
    plot(sim26b.xx, sim26b.vv * 3.6);
    hold off
    
    subplot(M,N,P); P = P + 1;
    plot(sim26a.xx, sim26a.stats.pbrakes * 1e-3);
    xlabel("Distance (m)")
    ylabel("Brakes (kW)")
    hold on
    plot(sim26b.xx, sim26b.stats.pbrakes * 1e-3);
    legend(["EV26" "EV26+"])
    hold off
    
    subplot(M,N,P); P = P + 1;
    plot(sim26a.xx, sim26a.stats.ptraction * 1e-3);
    xlabel("Distance (m)")
    ylabel("Tractive Power (kW)")
    hold on
    plot(sim26b.xx, sim26b.stats.ptraction * 1e-3);
    hold off
    
    subplot(M,N,P); P = P + 1;
    plot(sim26a.xx, sim26a.stats.pelectric * 1e-3);
    xlabel("Distance (m)")
    ylabel("Electric Power (kW)")
    hold on
    plot(sim26b.xx, sim26b.stats.pelectric * 1e-3);
    hold off
    
    subplot(M,N,P); P = P + 1;
    plot(sim26a.xx, sim26a.stats.ax / 9.806);
    xlabel("Distance (m)")
    ylabel("Longitudinal G Force (G)")
    ylim([-2 2])
    hold on
    plot(sim26b.xx, sim26b.stats.ax / 9.806);
    hold off
    
    subplot(M,N,P); P = P + 1;
    plot(sim26a.xx, sim26a.stats.ay / 9.806);
    xlabel("Distance (m)")
    ylabel("Lateral G Force (G)")
    ylim([-2 2])
    hold on
    plot(sim26b.xx, sim26b.stats.ay / 9.806);
    hold off
    figure(i+1)
    subplot(2,1,1);
    track_map(sim26a);
    title("EV26a")

    subplot(2,1,2); 
    track_map(sim26a);
    title("EV26b")
end

%% Simulation
function sim = summarize_lapsim(regs, track, car)
    sim.regs = regs;
    sim.track = track;
    sim.car = car;
    [sim.dx, sim.xx, sim.vv] = lapsim(sim.regs, sim.track, sim.car);
    [sim.dt, sim.tt, sim.splits] = laptime(sim.xx, sim.vv, sim.track);
    sim.time = sim.tt(end);
    sim.stats = calculate_power(sim.regs, sim.track, sim.car, sim.xx, sim.vv, sim.dt, sim.tt);
end

function [dx, xx, vv] = lapsim(regs, track, car)
    N = length(track.segments.lengths);
    cumlength = cumsum(track.segments.lengths);
    total_length = cumlength(end);
    dx = .005; % max distance to step by
    n = ceil(total_length / dx);
    dx = total_length / n; % ensure even division
    xx = 0:dx:total_length;
    gross_limits = arrayfun(@(k) calculate_limit(regs, track, car, k), track.segments.curvature);
    gross_limits = min(gross_limits, track.segments.limits); % by segments
    fine_limits = arrayfun(@(x) gross_limits(find_segment(cumlength, x)), xx); % by dx steps
    vv = fine_limits;
    % for i = 1:N
    %     if speed_limits(prev(i, N)) < speed_limits(i)
    %         x0 = cumlengths(prev(i, N));
    %         smear = evolve(+dx, cumlength, track, x0, speed_limits(prev(i, N)), vv);
    %         vv = min(vv, smear);
    %     elseif speed_limits(prev(i, N)) > speed_limits(i)
    %         smear = evolve(-dx, cumlength, track, x0, speed_limits(i), vv);
    %         vv = min(vv, smear);
    %     end
    % end
    for i = 1:n
        if fine_limits(prev(i, n)) < fine_limits(i)
            x0 = xx(i);
            v0 = fine_limits(prev(i, n));
            smear = evolve(+dx, cumlength, regs, track, car, i, x0, v0, vv);
            vv = min(vv, smear);
        elseif fine_limits(prev(i, n)) > fine_limits(i)
            x0 = xx(prev(i, n));
            v0 = fine_limits(i);
            smear = evolve(-dx, cumlength, regs, track, car, prev(i, n), x0, v0, vv);
            vv = min(vv, smear);
        end
    end
    % figure
    % plot(xx, fine_limits*3.6);
    % xlabel("Distance (m)")
    % ylabel("Velocity (km/h)")
    % hold on
    % plot(xx, vv * 3.6)
    % legend(["Max speed to complete sector" "Lap after time-evolution"])
end

function limit = calculate_limit(regs, track, car, k)
   %  top_speed = min((regs.power_limit / (0.5 * car.aero.cda * track.air_density))^(1/3), car.drive.motor.max_rpm * 2 * pi / 60 / car.drive.ratio * car.tires.radius);
   % 
   %  % 0.5 p cda v^2 = m a = m 0.5 mux (9.806 + .5 p cla v^2/m)
   %  % .5 p v^2 (cda - 0.5 mux cla) = 0.5 mux 9.806
   %  % v^2 = .5 mux 9.806 / (.5 p) / (cda - 0.5 mux cda)
   %  long_aero_limit = sqrt(9.806 * 0.5 * car.tires.mux / (0.5 * track.air_density) / (car.aero.cda - 0.5 * car.tires.mux * car.aero.cla));
   %  if imag(long_aero_limit) ~= 0
   %      long_aero_limit = Inf;
   %  end
   % 
   %  % v^2 k = muy (9.806 * mass + .5 p cla v^2) / mass
   %  % v^2 (k - .5 p cla muy / mass) = muy * g
   %  % v^2 = 9.806 / (|k|/muy - .5 p cla/mass)
   % 
   %  % todo: derate muy by fx from drag! otherwise this assumes no drag in
   %  % steady state cornering
   %  lat_aero_limit = sqrt(9.806 * car.tires.muy / (abs(k) - 0.5 * track.air_density * car.aero.cla * car.tires.muy / car.mass));
   %  if imag(lat_aero_limit) ~= 0
   %      lat_aero_limit = Inf;
   %  end
   %  lat_non_aero_limit = sqrt(9.806 * car.tires.muy / abs(k));
   % 
   %  options = [lat_aero_limit, top_speed];
   % % options = [lat_non_aero_limit, top_speed];
   %  limit = min(options);

    function accel = accel_from_vel(vel)
        [~, accel] = long_env(regs, track, car, k, v);
    end
    
    % turns out time-evolving is easier - just terminal velocity
    v = 0; a = 1; dt = 0.001;
    while abs(a) > 0.001 * 9.806
        a = accel_from_vel(v);
        v = v + a * dt;
    end
    limit = v;
end

function idx = find_segment(cumlength, x)
    idx = find(cumlength >= x, 1);
end

function p = prev(i, l)
    if i == 1
        p = l;
    else
        p = i - 1;
    end
end

function n = next(i, l)
    if i == l
        n = 1;
    else
        n = i + 1;
    end
end

function smear = evolve(dx, cumlength, regs, track, car, i0, x0, v0, limits)
    smear = limits;
    i = i0;
    di = sign(dx);
    v = v0;
    x = x0;
    while v <= limits(i)
        segment_idx = find_segment(cumlength, x);
        k = track.segments.curvature(segment_idx);
        [brake, accel] = long_env(regs, track, car, k, v);

        if sign(dx) == 1
            a = accel;
        else
            a = brake;
        end

        dv = a * dx / v;
        v = v + dv;

        smear(i) = v;

        if sign(dx) == 1
            i = next(i, length(limits));
        else
            i = prev(i, length(limits));
        end
    end
end

function [brake, accel] = long_env(regs, track, car, k, v)
    downforce = 0.5 * track.air_density * car.aero.cla * v^2;

    fz = 9.806 * car.mass + downforce;
    fy = v^2 * k * car.mass;
    effective_mux = sqrt(1 - min(1,(fy/(fz * car.tires.muy))^2)) * car.tires.mux;
    potential_fx = effective_mux * fz;

    drag = 0.5 * track.air_density * car.aero.cda * v^2;
    
    brake = (-potential_fx - drag) / car.mass;

    power_limited_accel = regs.power_limit / car.mass / v;
    % ax = (0.5 + h/l * ax/g) (potential_fx)/m
    % ax (1 - h/l potential/mg) = 0.5 potential / m
    % ax = 0.5 potential_fx / m / (1 - h/l potential_fx/mg)
    % ax = 0.5 / (m/potential - h/lg)
    traction_limited_accel = 0.5 * potential_fx / car.mass / (1 - car.cg(3) / 1000 / car.l * potential_fx/car.mass/9.806);
    
    w = v / car.tires.radius * car.drive.ratio;
    max_torque = interp1(car.drive.motor.LUT_WM_W, car.drive.motor.LUT_WM_M, w);
    motor_limited_accel = max_torque * car.drive.count * car.drive.ratio * car.drive.efficiency / car.tires.radius / car.mass;
    accel = min([power_limited_accel traction_limited_accel motor_limited_accel]) - drag / car.mass;
end

function [dt, tt, splits] = laptime(xx, vv, track)
    cumlength = cumsum(track.segments.lengths);    
    tt = 0 * xx;
    dt = 0 * xx;
    dx = xx(2) - xx(1); % constant
    splits = 0 * track.segments.lengths;
    for i = 1:length(xx)
        dt(i) = dx / vv(i);
        tt(i) = tt(prev(i, length(xx))) + dt(i);
        segment_idx = find_segment(cumlength, xx(i));
        splits(segment_idx) = tt(i);
    end
end

function stats = calculate_power(regs, track, car, xx, vv, dt, tt)
    cumlength = cumsum(track.segments.lengths);    
    stats.ptraction = xx * 0;
    stats.pbrakes = xx * 0;
    stats.pdrag = xx * 0;
    stats.pkinetic = xx * 0;
    stats.pelectric = xx * 0;
    stats.pmotor = xx * 0;  
    stats.tmotor = xx * 0;
    stats.wmotor = vv / car.tires.radius * car.drive.ratio;
    stats.ax = xx * 0;
    stats.ay = xx * 0;
    stats.k = xx * 0;

    for i = 1:length(xx)
        stats.ax(i) = (vv(next(i, length(vv))) - vv(i)) / (dt(i));
        segment_idx = find_segment(cumlength, xx(i));
        stats.k(i) = track.segments.curvature(segment_idx);
        stats.ay(i) = vv(i)^2*stats.k(i);
        f = car.mass * stats.ax(i);
        stats.pkinetic(i) = -f * vv(i);
        drag = 0.5 * car.aero.cda * track.air_density * vv(i)^2;
        stats.pdrag(i) = -drag * vv(i);
        ft = f + drag;
        if ft > 0
            stats.ptraction(i) = ft * vv(i);
            stats.pmotor(i) = stats.ptraction(i) / car.drive.efficiency;
            stats.tmotor(i) = ft * car.tires.radius / car.drive.ratio / car.drive.efficiency;
            stats.numotor(i) = interp2(car.drive.motor.LUT_WME_W, car.drive.motor.LUT_WME_M, car.drive.motor.LUT_WME_E.', clip(stats.wmotor(i), car.drive.motor.LUT_WME_W(1), car.drive.motor.LUT_WME_W(end)), clip(stats.tmotor(i), car.drive.motor.LUT_WME_M(1), car.drive.motor.LUT_WME_M(end)), 'nearest');
            stats.pelectric(i) = stats.pmotor(i) / stats.numotor(i);
        else
            stats.pbrakes(i) = ft * vv(i);
        end
    end

    stats.avgmotor = sum(dt .* stats.pmotor)/sum(dt);
    stats.avgelectric = sum(dt .* stats.pelectric)/sum(dt);
    stats.avgtraction = sum(dt .* stats.ptraction)/sum(dt);
    stats.avgbrakes = sum(dt .* stats.pbrakes)/sum(dt);
    stats.avgdrag = sum(dt .* stats.pdrag)/sum(dt);
end

function [xx, yy] = track_map(sim)
    x = 0;
    y = 0;
    theta = 0;

    xx = zeros(size(sim.vv));
    yy = zeros(size(sim.vv));
    
    for i=1:length(sim.vv)
        if sim.stats.k(i) == 0
            x = x + sim.dx * cos(theta);
            y = y + sim.dx * sin(theta);
        else
            radius = -1/sim.stats.k(i);
            delta_theta = -sim.dx * sim.stats.k(i);

            cx = x - radius * sin(theta);
            cy = y + radius * cos(theta);
            
            x = cx + radius * sin(theta + delta_theta);
            y = cy - radius * cos(theta + delta_theta);
            
            theta = theta + delta_theta;
        end

        xx(i) = x;
        yy(i) = y;
    end

    hold on;
    
    %I = imread('2025endurance_rot.png'); 
    %h = image([-29,51.5],[-4.5,6.5],I); 
    %uistack(h,'bottom')
    
    grid on;
    daspect([1 1 1])
    
    xticks(min(xx)-100:100:max(xx)+100);
    yticks(min(yy)-100:100:max(yy)+100);
    %fnplt(cscvn(xy), 'r', 2);
    
    % hold off;
    % axis equal;
    % 
    % %figure(2)
    % hold on;
    % axis equal;

    % for k=1:size(test_x,2)
    %     color = [cmap(k) 0 (1 - cmap(k))];
    %     curr_list_x = test_x{k};
    %     curr_list_y = test_y{k};
    %     plot(curr_list_x,curr_list_y,'.','Color',color);
    % end
    
    scatter(xx, yy, [], sim.vv.')
end

function [score] = lap2score(accel_time, skid_time, autox_time, endur_time, energy)
    %Best times from 2025 Michigan
    accel_min = 3.821;
    skidpad_min = 4.933;
    autocross_min = 45.734;
    endurance_min = 1369.936;

    score = 0;

    %Accel
    accel_max = 1.5*accel_min;
    if accel_time <= 0
        %
    elseif accel_time < accel_min
        score = score + 100;
    elseif accel_time < accel_max
        score = score + 4.5+95.5*((accel_max/accel_time)-1)/((accel_max/accel_min)-1);
    else
        score = score + 4.5;
    end

    %Skidpad
    skidpad_max = 1.25*skidpad_min;
    if skid_time <= 0
        %
    elseif skid_time < skidpad_min
        score = score + 75;
    elseif skid_time < skidpad_max
        score = score + 3.5 + 71.5*((skidpad_max/skid_time)^2-1)/((skidpad_max/skidpad_min)^2-1);
    else
        score = score + 3.5;
    end


    %Autocross
    autocross_max = 1.45*autocross_min;
    if autox_time <= 0
        %
    elseif autox_time < autocross_min
        score = score + 125;
    elseif autox_time < autocross_max
        score = score + 6.5+118.5*((autocross_max/autox_time)-1)/((autocross_max/autocross_min)-1);
    else
        score = score + 6.5;
    end

    %Endurance - Assumes finishing all 22 laps in endurance if a non-negative time
    endurance_max = 1.45*endurance_min;
    if endur_time <= 0
        %
    elseif endur_time < endurance_min
        score = score + 250+25;
    elseif endur_time < endurance_max
        score = score + 25+ 250*((endurance_max/endur_time)-1)/((endurance_max/endurance_min)-1);
    else
        score = score + 25;
    end
    
    e_factor_max = 0.848; %From UConn 2025
    e_factor_min = (62.270/(1.45*62.270))*((3.275*0.65/22)/(20.02*22/(100*22))); %0.483, Min time SJSU, Min CO2 UConn
    %Efficiency
    if endur_time <= 0
        
    elseif endur_time < endurance_min

    elseif endur_time < endurance_max
        e_factor = ((62.270)/(endur_time/22))*((3.275*0.65/22)/(energy/1000*0.65/22));
        e_score = 100*(e_factor-e_factor_min)/(e_factor_max-e_factor_min);
        e_score
        score = score + e_score;
    else
        
    end


end

