clear; clc; 

%% PARAMETERS
rho = 1240; c = 1800; k_th = 0.13;          % assuming material printed is PLA
rhoc  = rho*c;  kappa = k_th/rhoc;          % kappa' = k/(rho c) [m^2/s]
d = 0.4e-3;  A = d^2;                       % square cross-section [m], [m^2]
v = 40e-3;                                  % print speed [m/s]
a_ax = 20e-3; b_ax = 10e-3;                 % ellipse semi-axes [m]
L = pi*(3*(a_ax+b_ax) - sqrt((3*a_ax+b_ax)*(a_ax+3*b_ax)));  % approximate ellipse perimeter
T_noz = 220; T_plate = 60; T_amb = 25;      % [degC]
h_air = 10;  h_plate = 300;                 % [W/m^2K]

%% DERIVED QUANTITIES
a_air   = h_air*3*d/(A*rhoc);               % 3 exposed sides -> loss rate [1/s]
a_plate = h_plate*d/(A*rhoc);               % bottom side, on the plate [1/s]

dx = d/2;  N = ceil(L/dx);   % differential length size and number of steps for ellipse to be printed
dt = dx/v;                                  % exactly one new cell per time step
assert(kappa*dt/dx^2 <= 0.5, 'Diffusion (CFL) stability violated');
assert(dt*(a_air+a_plate) < 1, 'Time step too large for the side losses');

tol = 0.25;    % tolerance for nearly uniform temperature
t_max = (N + 10/(a_air+a_plate)/dt)*dt;     % Buffer T_max if ellipse not printed yet
K_max = ceil(t_max/dt);
n_save = 10;   % storing Temperature values every 20steps
n_store = ceil(K_max/n_save); % number of Temperature vectors stored

%% INITIALISATION
T = zeros(N,1);  n = 0;  % n = number of active cells
Tmax = zeros(1,K_max); Tmin = Tmax; Tmean = Tmax;
% Matrix to store temperature value across filament across time
T_hist = NaN(N,n_store); 
t_hist = zeros(1,n_store); % corresponding time for each temperature column
j = 0;

%% TIME MARCHING
for k = 1:K_max
    t = (k-1)*dt;

    % domain growth, one new cell added while printing (N)
    n_new = min(N,k);
    T(n+1:n_new) = T_noz;
    n = n_new;

    % Assuming Neumann conditions (no heat flux through filament ends).
    % Conduction occurs only for cells in between.
    q = zeros(n+1,1);                       
    q(2:n) = kappa*(T(1:n-1) - T(2:n))/dx;

    % Convective heat transfer between filament-plate and filament-air 
    f = a_air*(T_amb - T(1:n)) + a_plate*(T_plate - T(1:n));

    % Finite Difference Time Difference (FDTD) update 
    T(1:n) = T(1:n) + dt*((q(1:n) - q(2:n+1))/dx + f);

    Tmax(k) = max(T(1:n)); Tmin(k) = min(T(1:n)); Tmean(k) = mean(T(1:n));
    % saving the Temperature vector every 20 time steps
    if mod(k,n_save) == 0
        j = j+1; 
        T_hist(1:n,j) = T(1:n); % temperature update
        t_hist(j) = t; % correspoding time
    end
    
    % Time marching stopped if uniform temperature is achieved after
    % printing complete
    if k >= N && (Tmax(k) - Tmin(k)) < tol
        break; end
end
t_end = (k-1)*dt;
Tmax = Tmax(1:k); Tmin = Tmin(1:k); Tmean = Tmean(1:k);
T_hist = T_hist(:,1:j); t_hist = t_hist(1:j);
t_vec = (0:k-1)*dt;
t_print = (N-1)*dt;
fprintf('Printing ends at %.2f s, nearly uniform at %.2f s\n', t_print, t_end);

%% PLOTS
x_mm = ((1:N)-0.5)*dx*1e3;
figure; imagesc(t_hist, x_mm, T_hist, 'AlphaData', ~isnan(T_hist));
set(gca,'YDir','normal'); colorbar; xlabel('t [s]'); ylabel('path coordinate [mm]');
title('Temperature along the filament');

figure; plot(t_vec,Tmax,t_vec,Tmean,t_vec,Tmin); xline(t_print,'--');
xlabel('t [s]'); ylabel('T [^\circC]'); legend('max','mean','min','print ends');
