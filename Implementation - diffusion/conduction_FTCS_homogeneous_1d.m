%% 1D conduction - Forward time centered space
clc; clear; close all;
%% settings
% Parameters
L = 0.3;             % Length of the rod (m)
Nx = 300;             % Number of spatial points
x = linspace(0, L, Nx); % Discretized spatial domain

dx = x(2) - x(1);    % Spatial step size
time_max = 500;      % Total simulation time (s)
dt = 0.001;           % Time step (s)
Nt = time_max / dt; % Number of time steps

% Material properties (Aluminum)
k = 205;             % Thermal conductivity (W/m.K)
rho = 2700;          % Density (kg/m^3)
c = 897;             % Specific heat capacity (J/kg.K)
alpha = k / (rho * c); % Thermal diffusivity (m^2/s)

% Stability condition (Fourier number)
Fo = alpha * dt / dx^2;
if Fo > 0.5
    error('Unstable solution: Reduce dt or increase dx');
end

% Initial and Boundary Conditions
T = 25 * ones(Nx, 1); % Initial temperature (°C)
T(1) = 500;           % Left boundary condition (°C)
T(end) = 25;          % Right boundary condition (°C)


%% Time-stepping loop
for t = 1:Nt
    T_new = T; % temporary for new temperatures
    
    % Finite Difference Scheme (Explicit Method)
    for i = 2:Nx-1
        T_new(i) = T(i) + Fo * (T(i+1) - 2*T(i) + T(i-1));
    end
    
    T = T_new; % Update temperature distribution
    
    % Visualization - every 50 seconds
    if mod(t, 50000) == 0  
        plot(x, T, 'LineWidth', 2);
        hold on;
    end
end


%% plot
xlabel('Position along the rod (m)');
ylabel('Temperature (°C)');
title('Temperature Distribution Over Time');
grid on;
legend('t=50s', 't=100s', 't=150s', 't=200s', 't=250s', 't=300s', 't=350s', 't=400s', 't=450s', 't=500s');


