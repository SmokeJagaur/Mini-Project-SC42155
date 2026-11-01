clc; clear; close all;

%% Parameters
m_p = 0.5;  r = 0.01;  J = 5e-6;
b_theta = 1e-6;
r_p = 0.005;
c = b_theta/r_p^2;
b_linx = 0.2;
b_liny = 0.2;
beq_x = b_linx + 8*c;
beq_y = b_liny + 4*c;
M_eff = m_p + (2/r^2)*J;
printer = struct('r', r, 'Km', 0.05, 'meff', M_eff, 'bx', beq_x, 'by', beq_y,...
                 'Lm', 0.005, 'Rm', 5);   

%% Parameters for the square
L     = 0.2; % side length
P0    = [0.15; 0.15]; % starting point
T_seg = 3; 
dirs  = [1 0; 0 1; -1 0; 0 -1]';   % directions for each mode
% tolerances for the guards
tol_pos = 1e-6; % position
tol_vel = 1e-4; % velocity

%% Hybrid simulation
% z = [x y xd yd il ir thl_est thr_est]'
z  = [P0; 0; 0; 0; 0; 0; 0];  % starting with estimated angular positions as zero
t0 = 0;  P = P0;
opts = odeset('RelTol',1e-8, 'AbsTol',1e-10, 'MaxStep',T_seg/200);
T = [];  Z = [];  M = [];  t_sw = t0; U = [];

% iterating through modes, the next iteration occurs when the guards are
% detected
for k = 1:4                                   
    P1 = P + L*dirs(:,k); %end point of this line segment

    % system dynamics with input as the straight_line input function (using inverse dynamics)
    f = @(t,z) sys(z, stline_input(P(1),P(2),P1(1),P1(2),T_seg,t-t0,printer)', printer);

    % guards for switching modes
    opts.Events = @(t,z) guard(z, P0, P1, printer.r, tol_pos, tol_vel);

    % solving the mode's ODE using ODE45
    [t, zz] = ode45(f, [t0, t0 + 1.5*T_seg], z, opts);
    Uk = zeros(numel(t), 2);

    %storing the input values
    for j = 1:numel(t)
        Uk(j,:) = stline_input(P(1),P(2),P1(1),P1(2),T_seg,t(j)-t0,printer);
    end


    %logging outputs for plot
    T = [T; t];  Z = [Z; zz];  M = [M; k*ones(size(t))]; U = [U; Uk];
    z = zz(end,:)';  t0 = t(end);  P = P1;  t_sw(end+1) = t0;
end

% staying at the start after square is printed
[t, zz] = ode45(@(t,z) sys(z,[0;0],printer), [t0, t0+1], z, opts);
T = [T; t];  Z = [Z; zz];  M = [M; 5*ones(size(t))]; U = [U; zeros(numel(t),2)];


%% Plots
figure('Position',[100 100 900 750]);
subplot(2,1,1);
plot(T, Z(:,1), T, Z(:,2), 'LineWidth',1.5); hold on;
xline(t_sw(2:end), ':k');
grid on; xlabel('Time (s)'); ylabel('Position (m)'); legend('x','y');
title('Position');

subplot(2,1,2);
plot(T, Z(:,3), T, Z(:,4), 'LineWidth',1.5); hold on;
xline(t_sw(2:end), ':k');
grid on; xlabel('Time (s)'); ylabel('Velocity (m/s)'); legend('xd','yd');
title('Velocity');

figure('Name', 'Trajectory');
hold on;
names = {'q_r','q_u','q_l','q_d','q_s'};
for k = 1:5
    plot(Z(M==k,1), Z(M==k,2), 'LineWidth',1.8, 'DisplayName',names{k});
end
grid on; axis equal; xlabel('x (m)'); ylabel('y (m)');
xlim([0 0.4]); ylim([0 0.4]);
title('2D trajectory by mode'); legend('Location','eastoutside');

figure('Name','Motor voltages');
plot(T, U(:,1), T, U(:,2), 'LineWidth',1.5); hold on;
xline(t_sw(2:end), ':k');
xlim([0 19])
ylim([-1.2 1.2])
grid on; xlabel('Time (s)'); ylabel('Voltage (V)');
legend('u_l','u_r'); title('Motor input voltages');

%% functions
function dz = sys(z, u, p)
    % z = [x y xd yd il ir thl thr]',  u = [ul; ur]
    xd = z(3);  yd = z(4);  il = z(5);  ir = z(6);
    wl = (xd + yd)/p.r;                 
    wr = (xd - yd)/p.r;                 
    dz = [ xd;
           yd;
           (-p.bx*xd + (p.Km/p.r)*(il + ir))/p.meff;
           (-p.by*yd + (p.Km/p.r)*(il - ir))/p.meff;
           (u(1) - p.Rm*il - p.Km*wl)/p.Lm;     
           (u(2) - p.Rm*ir - p.Km*wr)/p.Lm;     
           wl;                                  
           wr ];
end

function [val, isterminal, direction] = guard(z, P0, P1, r, tol_pos, tol_vel)
    % guard is triggered when position and velocity are within the
    % tolerance
    % estimating position from the angular integrated angular positions
    p_est = P0 + (r/2)*[z(7) + z(8); z(7) - z(8)];

    % estimated printer head speed
    wl = (z(3) + z(4))/r;  wr = (z(3) - z(4))/r;
    speed = (r/2)*hypot(wl + wr, wl - wr);

    val = max(norm(p_est - P1) - tol_pos, speed - tol_vel);

    isterminal = 1;    % stop signal for ODE45
    direction  = -1;   % to make sure firing occurs when head approches
end