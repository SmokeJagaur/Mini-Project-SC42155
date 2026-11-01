clc; clear; close all;

%% Parameters
m_p = 0.5;  r = 0.01;  J = 2e-5;
M_eff = m_p + (2/r^2)*J;
printer = struct('r', r, 'Km', 0.05, 'meff', M_eff, 'bx', 1, 'by', 1, ...
                 'Lm', 0.005, 'Rm', 5, ...                
                 'Rn', 3.6, 'Ln', 1e-3, 'Cn', 1e-3, ...   
                 'Ct', 2, 'Rt', 5, 'Te', 25);             
T_print = 200;  % nozzle printing temperature            
%input to stay at that position
printer.u_hold = sqrt((T_print - printer.Te)*printer.Rn/printer.Rt);
printer.u_max  = 12;                         

%% Guards
g.tol_pos = 5e-4; % position tolerance
g.tol_vel = 2e-3; % velocity tolerance
g.T_wipe  = 2;    % time for excess filament to be extruded
g.T_load  = 3;    % time to load the new filament
g.T_low   = 100;  % nozzle temperature at which new filament is loaded
g.T_high  = 195;  % nozzle temperature to move back

%% Initial conditions
% Assuming the simulation starts as soon as the C of F command is give
P_start = [0.30; 0.25];   
v0      = [0.05; 0.03];  %printer head speed at printing position          
P_wipe  = [0; 0];                  
v_avg   = 0.10;                    
T_min   = 1;                       
names   = {'q_s','q_m','q_w','q_h','q_l','q_c','q_r'};
tmax_f  = [5 NaN 2*g.T_wipe 120 2*g.T_load 120 NaN];   % max timeouts per mode

% motor currents for inital state of zero acceleration
il0 = (r/(2*printer.Km))*(printer.bx*v0(1) + printer.by*v0(2));
ir0 = (r/(2*printer.Km))*(printer.bx*v0(1) - printer.by*v0(2));
uh  = printer.u_hold;
z   = [P_start; v0; il0; ir0;  T_print; uh/printer.Rn; uh;  0];

%% Hybrid simulation
t0 = 0;  
Psave = []; % saved position
t_sw = [];  T = [];  Z = [];  Q = []; U = [];
opts = odeset('RelTol',1e-6, 'AbsTol',1e-9, 'MaxStep',0.05);


% iterating through modes
for q = 1:7
    % points for straight line movement
    A = z(1:2);  B = z(1:2);
    if q == 2, B = P_wipe; end
    if q == 7, B = Psave;  end
    if q == 2 || q == 7
        Tr = max(norm(B-A)/v_avg, T_min);  tmax = 1.5*Tr; % duration of stline movement
    else
        Tr = 1;  tmax = tmax_f(q);
    end

    % ode for the mode
    f = @(t,z) sys(z, mode_input(q, t-t0, A, B, Tr, printer), printer);
    opts.Events = @(t,z) guard(z, q, B, g);

    % integrate until the guard fires (terminal event) or the timeout
    [t, zz] = ode45(f, [t0, t0 + tmax], z, opts);

    % logging inputs given through modes
    Uk = zeros(numel(t), 3); 
    for j = 1:numel(t)
        Uk(j,:) = mode_input(q, t(j)-t0, A, B, Tr, printer)';
    end
    
    T = [T; t];  Z = [Z; zz];  Q = [Q; q*ones(size(t))]; U = [U; Uk];
    z = zz(end,:)';  t0 = t(end);               

    % reset maps
    if q == 1, Psave = z(1:2); end          
    z(10) = 0;  % reseting timer
    t_sw(end+1) = t0;
end

%% Plots
figure('Name','Nozzle and inputs','Position',[150 150 900 750]);

% nozzle temperature state
subplot(3,1,1);
plot(T, Z(:,7), 'LineWidth',1.5); hold on; 
xlim([0 40]); ylim([75 210]);
yline(g.T_low,  '--b', 'T_{low}');
yline(g.T_high, '--r', 'T_{high}');
grid on; xlabel('Time (s)'); ylabel('T_n ({\circ}C)');
title('Nozzle temperature (dotted: mode switches)');

% input to the nozzle heater
subplot(3,1,2);
plot(T, U(:,3), 'LineWidth',1.5); 
xlim([0 37]); ylim([0 13]);
grid on; xlabel('Time (s)'); ylabel('u_n (V)');
title('Nozzle heater input voltage');

% input to the motors
subplot(3,1,3);
plot(T, U(:,1), T, U(:,2), 'LineWidth',1.5); 
xlim([0 37]); ylim([-2 2])
grid on; xlabel('Time (s)'); ylabel('Voltage (V)');
legend('u_l','u_r');
title('Motor input voltages');

%% functions
function dz = sys(z, u, p)
    % z = [x y xd yd il ir Tn in uc tm]',  u = [ul; ur; un]
    xd = z(3);  yd = z(4);  il = z(5);  ir = z(6);
    Tn = z(7);  in = z(8);  uc = z(9);
    wl = (xd + yd)/p.r;  % left motor angular velocity
    wr = (xd - yd)/p.r;  % right motor angular velocity
    dz = [ xd;
           yd;
           (-p.bx*xd + (p.Km/p.r)*(il + ir))/p.meff;
           (-p.by*yd + (p.Km/p.r)*(il - ir))/p.meff;
           (u(1) - p.Rm*il - p.Km*wl)/p.Lm;            
           (u(2) - p.Rm*ir - p.Km*wr)/p.Lm;            
           uc^2/(p.Rn*p.Ct) - (Tn - p.Te)/(p.Rt*p.Ct); 
           (u(3) - uc)/p.Ln;                           
           in/p.Cn - uc/(p.Rn*p.Cn);                   
           1 ];  % timer state
end

function u = mode_input(q, tl, A, B, Tr, p)
    % Input of mode q at current time
    uv = [0; 0]; % motor inputs for no motion modes
    if any(q == [2 7])                           
        uv = stline_input(A(1), A(2), B(1), B(2), Tr, tl, p)'; % stline input for move modes
    end
    un = p.u_hold;  % input for nozzle to hold temperature
    if q == 4 || q == 5
        un = 0;  % heater off when the filament is being loaded and turning off
    end             
    if q == 6     
        un = p.u_max; % heater on 
    end             
    u = [uv; un];
end

function [val, isterminal, direction] = guard(z, q, B, g)
    % Event function = guard of mode q. val is positive while the mode invariant
    % holds and crosses zero (downwards) when the guard becomes true.
    speed = sqrt(z(3)^2 + z(4)^2);
    % tolerance for st line moving modes
    B_reached  = max(norm(z(1:2) - B) - g.tol_pos, speed - g.tol_vel);
    if q == 1
        val = speed - g.tol_vel; % printer head stopped
    elseif q==2
        val = B_reached; % printer head reaches origin (wipe position)
    elseif q==3 
        val = g.T_wipe - z(10); % time for excess filament extrusion
    elseif q==4
        val = z(7) - g.T_low; % nozzle reaches filament loading temperature
    elseif q==5
        val = g.T_load - z(10); %loading time reached
    elseif q==6
        val = g.T_high - z(7); % printing temperature reached
    else
        val = B_reached; % printer head reaches saved position
    end

    isterminal = 1; % status for guard trigger
    direction  = -1;  % for ensuring trigger occurs when value reaches from a greater value
end