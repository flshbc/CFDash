%% FVM solution of 1D steady-state scalar transport equation
% flshbc @ Github, Spring 2025
% Associated function: ssst_eq_1d(.m)
%% Intro
% convective term = diffusive term 
% \int_S \rho \phi \bm{u} \cdot \bm{n} dS = \int_S \Gamma \nabla \phi \cdot \bm{n} dS 
% differential form: \nabla \cdot (\rho \phi \bm{u}) = \nabla \cdot (\Gamma \nabla \phi)
% unknown: \phi, others: constant, u -> Pe 
close all;clear;clc;
%% GE
rho = 1; % density of fluid
L = 1; % length of the region
gamma = 1; % convection coefficient
% u: constant velocity, used to control Pe
%% geo,BC,IC
% N: number of cells
phi_0 = 0;
phi_L = 1;
%% solving
u_test = 20;
N_test = 20;
[phi_linear,~,dx,xm] = ssst_eq_1d(rho,u_test,L,gamma,N_test,phi_0,phi_L,'linear-2nd','linear');
[phi_upwind,~,~] = ssst_eq_1d(rho,u_test,L,gamma,N_test,phi_0,phi_L,'upwind-2nd','upwind');

%% analytical solution
xm_demo = 0:0.001:1;
phi_exact = @(x,Pe,phi_0,phi_L,L) phi_0 + (exp(Pe.*x./L) - 1)./(exp(Pe) - 1)*(phi_L - phi_0);
%% plot solution
f1=figure();hold on;grid on; box on;
plot(xm_demo,phi_exact(xm_demo,u_test,phi_0,phi_L,L),'r-.','LineWidth',2);
plot(xm,phi_linear,'b+-','LineWidth',2);
plot(xm,phi_upwind,'gx-','LineWidth',2);
title(sprintf('1D steady-state scalar transport (N = %d, Pe = %.2f)',N_test,u_test));
xlabel('x'); ylabel('\phi');
legend('exact','numerical - linear','numerical - upwind');

%% error reuslts
eps_vector_linear = abs(phi_linear - phi_exact(xm,u_test,phi_0,phi_L,L)');
eps_linear = mean(eps_vector_linear);
eps_vector_upwind = abs(phi_upwind - phi_exact(xm,u_test,phi_0,phi_L,L)');
eps_upwind = mean(eps_vector_upwind);
fprintf('single case: N = %d, Pe = %.2f.\n',N_test,u_test);
fprintf('global relative error: %.3f%% (linear), %.3f%% (upwind).\n',100*eps_linear,100*eps_upwind);
%% order evaluation: (N - dx) power law
%  [eps ~（dx)^m]  equals to [log(eps) = m log(dx) + C]
u_order_test = 20;
% repeat solving and error computation for N - dx
N_set = 20:20:1000;
test_num = length(N_set);
dx_set = L./N_set;
eps_linear_loop = zeros(1,test_num);
eps_upwind_loop = zeros(1,test_num);
fprintf('\nmultiple case (%.d N): Pe = %.2f. \nestimated order of error: ',test_num,u_order_test);
for i = 1:test_num
    [phi_loop_linear,~,~,xm_loop] = ssst_eq_1d(rho,u_order_test,L,gamma,N_set(i),phi_0,phi_L,'linear-2nd','linear');
    [phi_loop_upwind,~,~,~] = ssst_eq_1d(rho,u_order_test,L,gamma,N_set(i),phi_0,phi_L,'upwind-2nd','upwind');
    exact_loop=phi_exact(xm_loop,u_order_test,phi_0,phi_L,L)';
    eps_linear_loop(i) = mean(abs(phi_loop_linear - exact_loop));
    eps_upwind_loop(i) = mean(abs(phi_loop_upwind - exact_loop));
end
% polyfit - linear fit (1st order) to get slope
p_linear = polyfit(log(dx_set),log(eps_linear_loop),1); 
p_upwind = polyfit(log(dx_set),log(eps_upwind_loop),1);
order_linear = p_linear(1);
order_upwind = p_upwind(1);
%% plot double-log
f2=figure();
loglog(dx_set,eps_linear_loop,'LineWidth',2);
hold on; grid on; box on;
loglog(dx_set,eps_upwind_loop,'LineWidth',2);
title('Error - grid size log-log relation');
xlabel('dx'); ylabel('log(\epsilon)');
legend('linear','upwind');
fprintf('%.4f (linear), %.4f (upwind).\n',order_linear,order_upwind);
