%% steady-state analytical solution for multi-layer material

total_length = 60;
n_layer = 2;
layer_choice = "uniform";

% layer
if layer_choice == "uniform"
    layer_length = ones(n_layer,1) * total_length / n_layer;
elseif layer_choice == "specified"
    layer_length = [];
    if sum(layer_length) ~= total_length
        error("incompatibale layer lengths");
    end
end

% 2 side dirichlet BC
t_left = 500;
t_right = 25;
% derive linear form 
t1 = @(x,c1,c2) (c1*x+c2);
t2 = @(x,c3,c4) (c3*x+c4);



k1 = 205; 
k2 = 30;
thick1 = 30;
thick2 = 5;

% balance heat flux crossing the layer q1 = q2
% q1 = k1*(T1 - T_inter)/thick1; % q = - k dt/dx;
% q2 = k2*(T_inter - T2)/thick2;
% q1 = q2;
q_uni = (t_left - t_right) / (thick1/k1 + thick2/k2);
% go from left to right, get T_inter(i)
T1 = t_left;
T2 = t_right;
T_inter_test = T1 - q_uni*thick1/k1;
T_inter = (k1/thick1*T1 +k2/thick2*T2)/(k1/thick1+k2/thick2);
% verify
q1 = k1*(T1-T_inter)/thick1;
q2 = k2*(T_inter - T2)/thick2;
if (abs(q1 - q2)<1e-5)
    q = q1;
end
c1 = -q/k1; % q = k*dt/dx ->  dt = q*dx/k   
c2 = t_left;
c3 = -q/k2;
c4 = T_inter- c3*thick1;

f2=figure();
hold on;
x1 = linspace(0,thick1,1000);
x2 = linspace(thick1,thick1+thick2,1000);
plot(x1,t1(x1,c1,c2),x2,t2(x2,c3,c4));
title('2-layer')