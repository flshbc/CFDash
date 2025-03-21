%% FVM solution of 1D unsteady scalar transport equation
% Associated function: usst_eq_1d(.m)
close all;clear;clc;
%% GE
% unsteady term = - convective term + diffusive term
rho = 1; % density of fluid
L = 1; % length of the region
gamma = 1; % diffusion coefficient
%% analytical solution
xm_demo = (0:0.001:1)';
phi_exact = @(x,Pe,phi_0,phi_L,L) phi_0 + ((exp(Pe.*x./L) - 1)./(exp(Pe) - 1)).*(phi_L - phi_0);
%% flow
u_test = 10; % constant velocity, used to control Pe
%Pe = rho*u_test*L/gamma;
%% spatial
N_test = 20; % number of cells
dx = L/N_test;
xm = linspace(0+dx/2,L-dx/2,N_test)'; % x-mid
%% temporal
% explained at function: usst_eq_1d
iter_max = 1000;
converge_crit = 1e-3;
stability_const = 1; % C
implicit_extent = 1; % theta
dt = stability_const/(abs(u_test)/dx + 2*gamma/(rho*dx^2));
CFL = abs(u_test)*dt/dx;
%% BC
phi_0 = 0; % dirichlet
phi_L = 1;
%% IC
dphi = (phi_L-phi_0)/N_test;
phi_ini = linspace(phi_0+dphi/2,phi_L-dphi/2,N_test); % assume no flow at the beginning
% phi_ini = zeros(N_test,1);
% phi_ini = phi_0*ones(N_test,1);
% phi_ini = phi_L*ones(N_test,1);
% phi_ini = - linspace(phi_0+dphi/2,phi_L-dphi/2,N_test);
% phi_ini = phi_exact(xm,u_test,phi_0,phi_L,L);

%% solving
[phi,A,Q,phi_diff,end_steps,end_time] = usst_eq_1d(rho,u_test,L,gamma,N_test, ...
    phi_0,phi_L,phi_ini, ...
    'linear-2nd','linear', ...
    iter_max,converge_crit,dt,implicit_extent);
%% plot
f1=figure('Position',[100,100,800,600],'Color','w');hold on;grid on;
time_slice = 5;
plot(xm_demo,phi_exact(xm_demo,u_test,phi_0,phi_L,L),'-','LineWidth',2,'DisplayName','Exact solution');
for k = 1:time_slice
    step_slice = round(k*end_steps/time_slice); % float to int
    plot([0; xm; L],[phi_0; phi(:,step_slice); phi_L], ...
        '-.','LineWidth',2, ...
        'DisplayName',sprintf('t = %.4f s',dt*(step_slice-1)));
end

title('FVM solution of 1D unsteady scalar transport');
xlabel('x');ylabel('\phi');
legend('Location','northwest');

%% usst_eq_1d: FVM solver for 1D unsteady scalar transport equation
% Associated file: unsteady_transport_FVM_1d.m
% reuse the steady-state solver for construction of const. A,Q
function [phi,A,Q,phi_diff,end_step,end_time] = usst_eq_1d(rho,u,L,gamma,N,phi_0,phi_L,phi_ini,bc,interp,iter_max,converge_crit,dt,theta)

% u: constant velocity, L: length of the region, gamma: convection coefficient, Pe = rho u L / gamma
% N: number of cells, iter_max: limit of total time steps
% phi, bc, interp: BCs, BC specs, interpolation specs

% stability_const (for dt): an overall strict time marching coefficient
%    based on Von Neumann analysis on explicit Euler method,
%    dt = const / (diff_lim + conv_lim) = const / (2*gamma/(rho*dx^2) + abs(u)/dx)
%    if const < 1, it ensures stability under both pure diff and conv
%    CFL is looser: Co = abs(u)*dt/dx

% theta: param in generalized temporal scheme theta-method, extent of algorithm being implicit
%    theta = 0   explicit Euler method
%    theta = 0.5 Crank-Nicolson
%    theta = 1   implicit Euler method

% converge_crit: the difference of solution between successive iterations is small enough
%    norm between old and new solutions is often used

% [1] preparation
dx = L/N;
A = zeros(N,N);
Q = zeros(N,1);
Pe = rho*u*L/gamma;
% [2] A phi = Q: boundary
if bc == "linear-1st"
    % linear/upwind is for interpolation of convective term at the boundary
    % 1st/2nd is for diffusive term (derivative) at the boundary
    A(1,1) = rho*u/2 + 3*gamma/dx;
    A(1,2) = rho*u/2 - gamma/dx;
    Q(1) = (rho*u + 2*gamma/dx)*phi_0;
    A(N,N) = -rho*u/2 + 3*gamma/dx;
    A(N,N-1) = -rho*u/2 - gamma/dx;
    Q(N) = (-rho*u + 2*gamma/dx)*phi_L;
elseif bc == "linear-2nd"
    A(1,1) = rho*u/2 + 4*gamma/dx;
    A(1,2) = rho*u/2 - 4/3*gamma/dx;
    Q(1) = (rho*u + 8/3*gamma/dx)*phi_0;
    A(N,N) = -rho*u/2 + 4*gamma/dx;
    A(N,N-1) = -rho*u/2 - 4/3*gamma/dx;
    Q(N) = (-rho*u + 8/3*gamma/dx)*phi_L;
elseif bc == "upwind-1st"
    A(1,1) = max(rho*u,0)-min(rho*u,0) + 3*gamma/dx;
    A(1,2) = min(rho*u,0) - gamma/dx;
    Q(1,1) = (max(rho*u,0) + 2*gamma/dx)*phi_0;
    A(N,N) = max(rho*u,0)-min(rho*u,0) + 3*gamma/dx;
    A(N,N-1) = -max(rho*u,0) - gamma/dx;
    Q(N) = (-min(rho*u,0) + 2*gamma/dx)*phi_L;
elseif bc == "upwind-2nd"
    A(1,1) = max(rho*u,0)-min(rho*u,0) + 4*gamma/dx;
    A(1,2) = min(rho*u,0) - 4/3*gamma/dx;
    Q(1,1) = (max(rho*u,0) + 8/3*gamma/dx)*phi_0;
    A(N,N) = max(rho*u,0)-min(rho*u,0) + 4*gamma/dx;
    A(N,N-1) = -max(rho*u,0) - 4/3*gamma/dx;
    Q(N) = (-min(rho*u,0) + 8/3*gamma/dx)*phi_L;
else
    error('Unexpected boundary condition setting (BC)');
end

% [3] A phi = Q: internal
if interp == "linear"
    for i = 2:N-1
        A(i,i-1) = -rho*u/2 - gamma/dx;
        A(i,i) = rho*u/2 - rho*u/2 + 2*gamma/dx;
        A(i,i+1) = rho*u/2 - gamma/dx;
    end
elseif interp == "upwind"
    for i = 2:N-1
        A(i,i-1) = -max([rho*u 0])-gamma/dx;
        A(i,i) = max([rho*u 0])-min([rho*u 0])+ 2*gamma/dx;
        A(i,i+1) = min([rho*u 0])-gamma/dx;
    end
else
    error('Unexpected interpolation scheme (interp)');
end

% [4] time marching: iteration
phi = zeros(N,iter_max);
phi_diff = zeros(iter_max,1);
phi(:,1) = phi_ini;
phi_diff(1) = Inf; % set large for first loop

LHS = eye(N,N)+theta*dt/dx*A;
RHS_matri = (eye(N,N)-(1-theta)*dt/dx*A);
RHS_Qterm = dt/dx*Q;

j = 1;
while(j <= iter_max && phi_diff(j) > converge_crit)
    RHS = RHS_matri*phi(:,j) + RHS_Qterm; % RHS contains previous info
    phi(:,j+1) = LHS \ RHS; % new A \ Q, A = I for explicit
    phi_diff(j+1) = norm(phi(:,j+1) - phi(:,j)); % absolute difference by norm
    j = j + 1;
end
end_step = j;
end_time = dt*(end_step-1);
fprintf("Pe = %.2f, N = %d, dt = %.2e  -->  iter %d, end time = %.2e s, diff = %.4e. \n",u,N,dt,end_step-1,end_time,phi_diff(end_step));
end

