%% FVM solution of 2D steady-state scalar transport equation
% flshbc @ Github | Spring 2025 | modification
%% Intro
% an extension based on fvm1dtransport.m, focused on specific 2D case
% inflow from south and west, while west BC is divided into 2 parts
% upwind/linear scheme for convection term, linear term for diffusion term
% study 1: artificial diffusion from upwind scheme (where gamma = 0 and diffusion still exists)
% study 2: linear scheme's unboundness/oscillatory behavior when |Pe_dx| > 2 (FVM solution violates maximum principle)
% study 3: indices and interpolations for 2D field
close all;clear;clc;
%% GE
rho = 1; % density
u = 1000; % velocity (forced convection), stands for Pe
v = 100;
Lx = 1;       % domain size
Ly = 1;
gamma = 1;  % diffusive coexfficient
%% scheme
conv='upwind'; %linear or upwind
%% geo
Nx = 40;      % Number of cells
Ny = 40;
Deltax = Lx/(Nx);
Deltay = Ly/(Ny);
N=Nx*Ny;

eq=@(i,j) j+(i-1)*Ny; % index convert
A = zeros(N,N);
b   = zeros(N,1);
%% BC
% BC setting
phi_s=0;   % south and west: Dirichlet, north and east: zero gradient
phi_w1=0;
phi_w2=1; 
% BC assigning
phi_w=zeros(Ny,1);
phi_w(1:int64(Ny/4),1)=phi_w1;
phi_w(int64(Ny/4)+1:end,1)=phi_w2;
%% dimensionless number
% Peclet(dx) numbers
Pe = rho*u*Lx/gamma;
Pe_y = rho*v*Ly/gamma;
Pe_dx=rho*u*Deltax/gamma;
Pe_dy=rho*v*Deltay/gamma;
fprintf('Nx = %d, Ny = %d\n',Nx,Ny);
fprintf('memory usage: %.4f MB\n',N*N/10^6*8);
fprintf('Pe_dx = %.4f, Pe_dy = %.4f\n',Pe_dx,Pe_dy);
%% FVM: upwind & linear 
switch conv
    case ('upwind')
        disp('convective term: upwind scheme');
        for i = 1:Nx
            for j = 1:Ny
                % West flux
                % Boundary cell
                if i==1
                    dif=2*gamma/Deltax*Deltay;
                    con=rho*u*Deltay;
                    b(eq(i,j),1) =  (dif +con)*phi_w(j);
                    A(eq(i,j),eq(i,j)) = A(eq(i,j),eq(i,j)) +dif;
                    % Internal cells
                else
                    dif=gamma/Deltax*Deltay;
                    con=rho*u*Deltay;
                    A(eq(i,j),eq(i,j)-Ny) = A(eq(i,j),eq(i,j)-Ny) -dif -max(con,0);
                    A(eq(i,j),eq(i,j)) = A(eq(i,j),eq(i,j)) +dif -min(con,0);
                end

                % East flux
                % internal
                if i<Nx
                    dif=gamma/Deltax*Deltay;
                    con=rho*u*Deltay;
                    A(eq(i,j),eq(i,j)+Ny) = A(eq(i,j),eq(i,j)+Ny) -dif + min(con,0);
                    A(eq(i,j),eq(i,j)) = A(eq(i,j),eq(i,j)) +dif +max(con,0);
                    % boundary
                elseif i==Nx
                    con=rho*u*Deltay;
                    A(eq(i,j),eq(i,j)) = A(eq(i,j),eq(i,j)) +max(con,0);
                end

                % South flux
                % boundary
                if j==1
                    dif=2*gamma/Deltay*Deltax;
                    con=rho*v*Deltax;
                    b(eq(i,j),1) =  (dif +con)*phi_s;
                    A(eq(i,j),eq(i,j)) = A(eq(i,j),eq(i,j)) +dif;
                    % internal
                else
                    dif=gamma/Deltay*Deltax;
                    con=rho*v*Deltax;
                    A(eq(i,j),eq(i,j)-1) = A(eq(i,j),eq(i,j)-1) -dif -max(con,0);
                    A(eq(i,j),eq(i,j)) = A(eq(i,j),eq(i,j)) +dif - min(con,0);
                end

                % North flux
                % internal
                if j<Ny
                    dif=gamma/Deltay*Deltax;
                    con=rho*v*Deltax;
                    A(eq(i,j),eq(i,j)+1) = A(eq(i,j),eq(i,j)+1) -dif + min(con,0);
                    A(eq(i,j),eq(i,j)) = A(eq(i,j),eq(i,j)) +dif + max(con,0);
                    %boundary
                elseif j==Ny
                    con=rho*v*Deltax;
                    A(eq(i,j),eq(i,j)) = A(eq(i,j),eq(i,j)) +con;
                end
            end
        end

    case {'linear'}
        disp('convective term: linear scheme');
        for i = 1:(Nx)
            for j = 1:Ny
                % West flux
                % Boundary
                if i==1
                    dif=2*gamma/Deltax*Deltay;
                    con=rho*u*Deltay;
                    b(eq(i,j),1) =  (dif +con)*phi_w(j);
                    A(eq(i,j),eq(i,j)) = A(eq(i,j),eq(i,j)) +dif;
                    % Internal
                else
                    dif=gamma/Deltax*Deltay;
                    con=rho*u*Deltay/2;
                    A(eq(i,j),eq(i,j)-Ny) = A(eq(i,j),eq(i,j)-Ny) -dif -con;
                    A(eq(i,j),eq(i,j)) = A(eq(i,j),eq(i,j)) +dif -con;
                end

                % East flux
                % internal
                if i<Nx
                    dif=gamma/Deltax*Deltay;
                    con=rho*u*Deltay/2;
                    A(eq(i,j),eq(i,j)+Ny) = A(eq(i,j),eq(i,j)+Ny) -dif +con;
                    A(eq(i,j),eq(i,j)) = A(eq(i,j),eq(i,j)) +dif +con;
                    % boundary
                elseif i==Nx
                    % dif=2*gamma/Deltax*Deltay; not used for outflow BC
                    con=rho*u*Deltay;
                    A(eq(i,j),eq(i,j)) = A(eq(i,j),eq(i,j)) +con;
                end

                % South flux
                % boundary
                if j==1
                    dif=2*gamma/Deltay*Deltax;
                    con=rho*v*Deltax;
                    b(eq(i,j),1) =  (dif +con)*phi_s;
                    A(eq(i,j),eq(i,j)) = A(eq(i,j),eq(i,j)) +dif;
                    % internal
                else
                    dif=gamma/Deltay*Deltax;
                    con=rho*v*Deltax/2;
                    A(eq(i,j),eq(i,j)-1) = A(eq(i,j),eq(i,j)-1) -dif -con;
                    A(eq(i,j),eq(i,j)) = A(eq(i,j),eq(i,j)) +dif - con;
                end

                % North flux
                % internal
                if j<Ny
                    dif=gamma/Deltay*Deltax;
                    con=rho*v*Deltax/2;
                    A(eq(i,j),eq(i,j)+1) = A(eq(i,j),eq(i,j)+1) -dif +con;
                    A(eq(i,j),eq(i,j)) = A(eq(i,j),eq(i,j)) +dif +con;
                    % boundary
                elseif j==Ny
                    dif=2*gamma/Deltay*Deltax;
                    con=rho*v*Deltax;
                    A(eq(i,j),eq(i,j)) = A(eq(i,j),eq(i,j)) +con;
                end
            end
        end

end  

phi = A\b; %[A][phi]=[b], dimension: N=Nx*Ny

%% data processing

x=[0 Deltax/2:Deltax:Lx-Deltax/2 Lx]; 
y=[0 Deltay/2:Deltay:Ly-Deltay/2 Ly];
phi_2D = zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        phi_2D(i,j)=phi(eq(i,j));
    end
end
phi_plot=zeros(Nx+2,Ny+2);             % result for plotting (incl. cell centre and boundary face values: 2 extra elements in each direction)
phi_plot(2:end-1,2:end-1)=phi_2D(:,:); % replace cell centre values from solution
phi_plot(:,1)=phi_s;                   % south and west boundary: Dirichlet
phi_plot(1,2:end-1)=phi_w(:);
phi_plot(:,end)=phi_plot(:,end-1);     % east and north boundaries: zero gradient
phi_plot(end,:)=phi_plot(end-1,:);

%% plot
f1 = figure('Units','normalized','Position',[0.1 0.25 0.8 0.5]);
subplot(1,2,1);
contourf(x,y,phi_plot',20,'edgecolor','none');
colormap('default');
hold on;
plot([0.5 0.5],[0 1],'r--');
set(gca,'Fontsize',14);
xlabel('x');
ylabel('y');
title('FVM solution');
axis equal; colorbar;
%clim([0,1]); % filter possible unboundness
subplot(1,2,2);
plot(y,phi_plot(Nx/2,:),'ro-','linewidth',2);
set(gca,'Fontsize',14);
xlabel('y');
ylabel('\phi');
title('solution on the mid-line');
hold on;grid on;axis equal;

