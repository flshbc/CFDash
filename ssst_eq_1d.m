%% FVM solver for 1D steady-state scalar transport equation
% flshbc, Spring 2025
% Associated file: fvm1dtransport.m
%% declaration
function [phi,Pe,dx,xm] = ssst_eq_1d(rho,u,L,gamma,N,phi_0,phi_L,bc,interp)
% u: constant velocity
% L: length of the region
% gamma: convection coefficient
% Pe = rho u L / gamma
% N: number of cells
%% parameters
Pe = rho*u*L/gamma; % Pe = rho u L / gamma
dx = L/N;
xm = linspace(0+dx/2,L-dx/2,N);
A = zeros(N,N);
phi = zeros(N,1);
Q = zeros(N,1);
%% A phi = Q: boundary
% FD/BD for 1 point/2 points
if bc == "linear-1st" % linear interpolation for advective term, 1st order FD/BD for diffusive term
    % diff: 1st-order FD/BD + conv: linear interp
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
    error('Unexpected boundary condition (BC) setting');
end
%% A phi = Q: interval
% 
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
%% solve
phi = A \ Q; % phi = A^{-1} Q
end