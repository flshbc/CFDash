% MATLAB script to animate the 1D Heat Equation solution (Heat Kernel)
clear; clc; close all;

%% Parameters
A0 = 1;                 % Amplitude of initial heat distribution
alpha = 1;             % Thermal diffusivity
L = 1;                 % Initial spread width (acts as delta approx)
k = 2 * alpha;         % Diffusion coefficient
x = linspace(-5, 5, 200);  % Space domain
t_max = 5;             % Maximum simulation time
dt = 0.05;             % Time step
t_vals = 0:dt:t_max;   % Time array

%% Set Up Figure for Animation
figure;
hold on;
h = plot(x, zeros(size(x)), 'b', 'LineWidth', 2); % Initialize plot
xlabel('x');
ylabel('Temperature u(x,t)');
title('Heat Kernel Solution Animation');
grid on;
ylim([0, 0.5]);  % Fixed y-axis limits for better visualization
xlim([-5, 5]);   % Fixed x-axis limits
legend_text = text(-4, 0.45, '', 'FontSize', 12, 'FontWeight', 'bold');

%% Animation Loop
for t = t_vals
    % Compute the heat kernel solution
    u_xt = (A0 / sqrt(2 * pi * (L^2 + k*t))) * exp(-x.^2 / (2 * (L^2 + k*t)));
    
    % Update plot
    set(h, 'YData', u_xt); 
    set(legend_text, 'String', sprintf('t = %.2f', t)); % Update time label
    
    % Pause for animation effect
    pause(0.05);
end

hold off;
