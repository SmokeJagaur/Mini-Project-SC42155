%% Authors 
% Sourish, Daniel

clc;
clear all;

%% Variables 

%printer head
m_p = 0.5; % mass printer head
%J_left = 0; % left motor disk inertia
%J_right = J_left; % right motor disk inertia
r = 0.01; % radius of motor disk

J = 2e-5; % assuming both motor disks to have same inertia
M_eff = m_p + (2/r^2)*J;
b_x = 1; %linear bearing friction x-direction
b_y = 1; %linear bearing friction y-direction
K_m = 0.05; %transducer electrical-mechanical

% motor circuit
R_m = 2000; % resistance
L_m = 0.005; % Inductance



%% state space 

A = [[0 0 1 0 0 0]
    [0 0 0 1 0 0]
    [0 0 -(b_x/M_eff) 0 K_m/(M_eff*r) K_m/(M_eff*r)]
    [0 0 0 -(b_y/M_eff) K_m/(M_eff*r) -K_m/(M_eff*r)]
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

%% Input vectors
% Input for only x direction motion
%u = [5*ones(size(t)); 5*ones(size(t))]'; 

% Input for only y direction motion
%u = [5*ones(size(t)); -5*ones(size(t))]'; 

u = zeros(length(t), 2);

% input as function of time
u1 = ((1/3)*sin(t))';
u2 = ((2/3)*cos(t))';

% final inputs
u = [u1, u2];


%% Simulating


[y, t_out, x] = lsim(sys, u', t);

%% Plotting Results
figure('Name', 'Printer Head Simulation', 'Position', [100, 100, 900, 600]);

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
