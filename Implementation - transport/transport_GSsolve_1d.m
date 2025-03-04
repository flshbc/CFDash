%% Gauss-Seidel method test on FVM solution of 1D steady-state scalar transport equation
% Associated function: ssst_eq_1d(.m) GSsolve(.m)
%% Intro
% convective term = diffusive term 
% \int_S \rho \phi \bm{u} \cdot \bm{n} dS = \int_S \Gamma \nabla \phi \cdot \bm{n} dS 
% differential form: \nabla \cdot (\rho \phi \bm{u}) = \nabla \cdot (\Gamma \nabla \phi)
% \phi is unknown, others are constant, u stands for Pe.

% Implement Gauss-Seidel method to solve the linear system, and compare the
% convergence between different relaxation factors.
close all;clear;clc;
%% GE
rho = 1; % density of fluid
L = 1; % length of the region
gamma = 1; % convection coefficient
% u: constant velocity, used to control Pe
% N: number of cells
%% BC
phi_0 = 0;
phi_L = 1;
%% construting & solving - matlab
u_test = 20;
N_test = 100;
fprintf('single case: N = %d, Pe = %.2f.\n',N_test,u_test);
[phi_linear,Pe,dx,xm,A,Q] = ssst_eq_1d(rho,u_test,L,gamma,N_test,phi_0,phi_L,'linear-2nd','linear');
%[phi_upwind,~,~,~,A2,Q2] = ssst_eq_1d(rho,u_test,L,gamma,N_test,phi_0,phi_L,'upwind-2nd','upwind');
%% analytical solution
xm_demo = 0:0.001:1;
phi_exact = @(x,Pe,phi_0,phi_L,L) phi_0 + (exp(Pe.*x./L) - 1)./(exp(Pe) - 1)*(phi_L - phi_0);
exact_phi = phi_exact(xm,Pe,phi_0,phi_L,L);
%% solving - Gauss Seidel method with relaxation
itermax=5000;
conv_abs_res = 1e-5;
[phi_gs,iter,res,err] = GSsolve(N_test,A,Q,itermax,conv_abs_res,1.0,exact_phi);
[phi_gs_o,iter_o,res_o,err_o] = GSsolve(N_test,A,Q,itermax,conv_abs_res,1.5,exact_phi); % over relax
[phi_gs_u,iter_u,res_u,err_u] = GSsolve(N_test,A,Q,itermax,conv_abs_res,0.5,exact_phi); % under relax

%% plot solution
f1=figure();hold on;grid on;box on;
plot(xm_demo,phi_exact(xm_demo,Pe,phi_0,phi_L,L),'r-.','LineWidth',2);
plot(xm,phi_linear,'b+-','LineWidth',2);
plot(xm,phi_gs,'gx-','LineWidth',2);
title(sprintf('1D steady-state scalar transport (N = %d, Pe = %.2f)',N_test,u_test));
xlabel('x'); ylabel('\phi');
legend('exact','matlab','GS');

%% error reuslts
eps_vector_linear = abs(phi_linear - exact_phi');
eps_linear = mean(eps_vector_linear);
%eps_vector_upwind = abs(phi_upwind - exact_phi');
%eps_upwind = mean(eps_vector_upwind);
eps_vector_gs = abs(phi_gs - exact_phi');
eps_gs = mean(eps_vector_gs);
fprintf('Finite volume method:\nglobal relative error: %.3f%% (matlab mldivide), %.3f%% (Gauss-Seidel, w=1).\n',100*eps_linear,100*eps_gs);

%% plot error and resiual vs iter
f2=figure();
loglog(1:iter_u,err_u(1:iter_u),'LineWidth',2);
hold on;grid on;box on;
loglog(1:iter,err(1:iter),'LineWidth',2);
loglog(1:iter_o,err_o(1:iter_o),'LineWidth',2);
title('Error vs num iter - Gauss-Seidel method');
xlabel('number of iteration'); ylabel('error');
legend('under relaxtion','\default','over relaxation');

f3=figure();
loglog(1:iter_u,res_u(1:iter_u),'LineWidth',2);
hold on;grid on;box on;
loglog(1:iter,res(1:iter),'LineWidth',2);
loglog(1:iter_o,res_o(1:iter_o),'LineWidth',2);
title('Residual vs num iter - Gauss-Seidel method');
xlabel('number of iteration'); ylabel('residual');
legend('under relaxtion','\default','over relaxation');

%% ssst_eq_1d: FVM solver(construction) for 1D steady-state scalar transport equation
% Associated file: transport_GSsolve_1d.m
%% declaration
function [phi,Pe,dx,xm,A,Q] = ssst_eq_1d(rho,u,L,gamma,N,phi_0,phi_L,bc,interp)
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

%% GSsolve: Gauss-Seidel method with relaxation
% Associated file: transport_GSsolve_1d.m
%% declaration
function [phi_gs,iter,res,err] = GSsolve(N,A,Q,itermax,conv_abs_res,w,exact_sol)
% N: Dimension of unknown phi
% A: N*N Q:N*1 (A*phi = Q)
% w: relaxation factor
%% preparations
phi_gs = zeros(N,1); % solution
res = zeros(itermax,1);
res(1) = norm(Q-A*phi_gs,Inf);
err = zeros(itermax,1);
err(1) = norm(exact_sol'-phi_gs,Inf);
iter = 1; % iteration count (1 shift)
%% iteration
while((res(iter) > conv_abs_res && iter < itermax))
    for i = 1: N
        I = [1:i-1 i+1:N];% two parts of indices to be operated on: [1:i-1 i+1:N]
        phi_gs(i) = w*(Q(i)-dot(A(i,I),phi_gs(I)))/A(i,i) + (1-w)*phi_gs(i); 
    end
    iter = iter + 1;
    res(iter) = norm(Q - A * phi_gs,Inf); % generally: Inf-norm < 2-norm < 1-norm
    err(iter) = norm(phi_gs - exact_sol',Inf);
end
%% result
fprintf('Gauss-Seidel method: \nrelaxation factor = %.2f, iteration = %d, absolute residual = %.3e\n',w,iter-1,res(iter));  
end