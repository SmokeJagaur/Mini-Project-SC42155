%% Authors 
% Sourish, Daniel

clc; clear; close all;

%% Variables 

%printer head
m_p = 0.5; % mass printer head
%J_left = 0; % left motor disk inertia
%J_right = J_left; % right motor disk inertia
r = 0.01; % radius of motor disk
r_p = 0.005; % idle pulley radius 
b_theta = 1e-6; % viscous friction for the bearing
c = b_theta/r_p^2; %equivalent belt damping for each pulley bearing

% motor and bearings
J = 5e-6; % assuming both motor disks to have same inertia
M_eff = m_p + (2/r^2)*J;

b_linx = 0.2; %linear bearing friction x-direction
b_liny = 0.2; %linear bearing friction y-direction

beq_x = b_linx + 8*c; % all bearings rotate for movement along x
beq_y = b_liny + 4*c; % bearings on the printer gantry do not rotate

K_m = 0.05; %transducer electrical-mechanical

% motor circuit
R_m = 5; % resistance
L_m = 5e-3; % Inductance

% in structure form
printer = struct('r', r, 'Km', K_m, 'meff', M_eff, 'bx', beq_x, 'by', beq_y, ...
           'Lm', L_m, 'Rm', R_m);
%% Constraints
dx = 0.5;
dy = 0.5;

%% state space 

A = [[0 0 1 0 0 0]
    [0 0 0 1 0 0]
    [0 0 -(beq_x/M_eff) 0 K_m/(M_eff*r) K_m/(M_eff*r)]
    [0 0 0 -(beq_y/M_eff) K_m/(M_eff*r) -K_m/(M_eff*r)]
    [0 0 -K_m/(r*L_m) -K_m/(r*L_m) -R_m/L_m 0]
    [0 0 -K_m/(r*L_m) K_m/(r*L_m) 0 -R_m/L_m]];

B = [[0 0]
    [0 0]
    [0 0]
    [0 0]
    [1/L_m 0]
    [0 1/L_m]];

C = eye(6);
D = zeros(6,2);

sys = ss(A, B, C, D);

t = 0:0.01:20; % time bounds
steps = size(t, 2);

%% Input vectors
%% straight line trajectory
% x1 = 0; y1 = 0;
% x2 = 0.25; y2 = 0.25;
% T = 3; % movement duration
%time
% t = 0:0.001:(T+1);
% steps = numel(t);

%% ellipse trajectory input
a = 0.2; b = 0.1; %semi-major and semi-minor axis of the ellipse
% drawing of the ellipse is split into three parts
% (1) origin to center straight line
T1 = 4;
x1 = 0; y1 = 0;
% (2) center to major axis end
T2 = 3;
x2 = dx/2; y2 = dy/2;
% (3) printing of the ellipse
T3 = 12;
x3 = x2 + a; y3 = y2;

% concatenated input
T = T1 + T2 + T3 + 1;
t = 0:0.001:T;
steps = numel(t);

u = zeros(steps, 2);
xr = zeros(steps, 1);  yr = zeros(steps, 1);   % reference path

for k = 1:steps
    tk = t(k);
    if tk < T1
        [u(k,:), ref] = stline_input(x1, y1, x2, y2, T1, tk, printer);
        xr(k) = ref.x;  yr(k) = ref.y;
    elseif tk < T1 + T2
        [u(k,:), ref] = stline_input(x2, y2, x3, y3, T2, tk - T1, printer);
        xr(k) = ref.x;  yr(k) = ref.y;
    else
        [u(k,:), ref] = ellipse_input(tk - T1 - T2, a, b, T3, printer);
        xr(k) = x2 + ref.x;  yr(k) = y2 + ref.y;   % shift ellipse to center
    end
end

%simulation
x0 = [x1; y1; 0; 0; 0; 0];
[y, t_out] = lsim(sys, u, t, x0);

%% Plotting Results
figure('Name', 'Printer Head Simulation', 'Position', [100, 100, 600, 600]);

% Subplot 1: Positions v/s time
subplot(3,1,1);
plot(t_out, y(:,1), 'LineWidth', 1.5, 'DisplayName', 'X Position');
hold on;
plot(t_out, y(:,2), 'LineWidth', 1.5, 'DisplayName', 'Y Position');
grid on;
xlabel('Time (s)');
ylabel('Position (m)');
title('Position vs Time for X and Y directions');
legend('Location', 'best');

%Subplot 2: Velocity v/s Time
subplot(3,1,2);
plot(t_out, y(:,3), 'LineWidth', 1.5, 'DisplayName', 'X velocity');
hold on;
plot(t_out, y(:,4), 'LineWidth', 1.5, 'DisplayName', 'Y velocity');
grid on;
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity vs Time for X and Y directions');
legend('Location', 'best');

% Subplot 3: 2D Trajectory (X vs Y)
subplot(3,1,3);
plot(y(:,1), y(:,2), 'Color', 'b', 'LineWidth', 1.5);
grid on;
xlabel('X Position (m)');
ylabel('Y Position (m)');
title('2D Trajectory of Printer Head');
axis equal;

% plot for motor input voltage and current
figure('Name', 'Motor Inputs and Currents', 'Position', [150 150 600 500]);

subplot(2,1,1);
plot(t, u(:,1), 'LineWidth', 1.5, 'DisplayName', 'u_{left}'); hold on;
plot(t, u(:,2), 'LineWidth', 1.5, 'DisplayName', 'u_{right}');
grid on; xlabel('Time (s)'); ylabel('Voltage (V)');
title('Motor Voltages'); legend('Location', 'best');

subplot(2,1,2);
plot(t_out, y(:,5), 'LineWidth', 1.5, 'DisplayName', 'i_{left}'); hold on;
plot(t_out, y(:,6), 'LineWidth', 1.5, 'DisplayName', 'i_{right}');
grid on; xlabel('Time (s)'); ylabel('Current (A)');
title('Motor Currents'); legend('Location', 'best');

