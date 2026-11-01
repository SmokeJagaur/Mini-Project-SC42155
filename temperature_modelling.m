clear; clc; 

%% parameters
rho = 1240; c = 1800; k_th = 0.13; % assuming material printed is PLA
rhoc  = rho*c;  kappa = k_th/rhoc;          
d = 0.4e-3;  A = d^2;                       
v = 40e-3;                                 
a_ax = 20e-3; b_ax = 10e-3;                
L = pi*sqrt(2*(a_ax^2 + b_ax^2));  % approximate ellipse perimeter
T_noz = 220; T_plate = 60; T_amb = 25;  
h_air = 10;  h_plate = 300;           

a_air   = h_air*3*d/(A*rhoc);         
a_plate = h_plate*d/(A*rhoc);            

dx = d/2;  N = ceil(L/dx);   % length step across the filament ellipse
dt = dx/v;  % time-step calculation based on speed of printing
assert(kappa*dt/dx^2 <= 0.5, 'Diffusion (CFL) stability violated');

tol = 0.5;    % tolerance for nearly uniform temperature
t_max = (N + 10/(a_air+a_plate)/dt)*dt;     % Buffer T_max if ellipse not printed yet
K_max = ceil(t_max/dt);
n_save = 20;   % storing Temperature values every 20steps
n_store = ceil(K_max/n_save); % number of Temperature vectors stored

%% initialization
T = zeros(N,1);  n = 0;  % n = number of active cells
Tmax = zeros(1,K_max); Tmin = Tmax; Tmean = Tmax;
% Matrix to store temperature value across filament across time
T_hist = NaN(N,n_store); 
t_hist = zeros(1,n_store); % corresponding time for each temperature column
j = 0;

%% Time marchine
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

figure('Name', 'Temperature Heat Map T(s, t)', 'Position', [500, 500, 600, 500]); 
imagesc(t_hist, x_mm, T_hist, 'AlphaData', ~isnan(T_hist));
set(gca,'YDir','normal'); colorbar; xlabel('t [s]'); ylabel('path coordinate [mm]');
title('Temperature along the filament');

figure('Name', 'Temperature Estimates vs Time', 'Position', [500, 500, 600, 500]); 
title('Temperature Eestimates vs Time')
plot(t_vec,Tmax,t_vec,Tmean,t_vec,Tmin, LineWidth=1.5); xline(t_print,'--');
xlim([0 17])
ylim([50 225])
xlabel('t [s]'); ylabel('T [^\circC]'); 
legend('max','mean','min','print ends');
