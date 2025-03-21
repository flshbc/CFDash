%% plot
f1=figure('Position',[100,100,800,600],'Color','w');hold on;grid on; box off;
time_slice = 5;
plot(xm_demo,phi_exact(xm_demo,u_test,phi_0,phi_L,L),'-','LineWidth',2,'DisplayName','Exact solution');
for k = 1:time_slice
    step_slice = round(k*end_steps/time_slice); % float to int
    plot([0; xm; L],[phi_0; phi(:,step_slice); phi_L], ...
        '-.','LineWidth',2, ...
        'DisplayName',sprintf('t = %.4f s',dt*(step_slice-1)));
end

ax=gca;
ax.FontName='default';ax.FontSize=12;
ax.FontWeight="bold";ax.FontAngle="italic";

title(sprintf('FVM solution of 1D unsteady scalar transport'),'FontAngle','normal');
xlabel('x','FontSize',13.2);ylabel('\boldmath$\phi$','Interpreter','latex','FontSize',13.2);
legend('Location','northwest');
xticks(0:L/5:L);yticks(phi_0:(phi_L-phi_0)/5:phi_L);

%exportgraphics(f1,'a.png','Resolution',300);
%cla;clf;
%% animation (standard)
f2=figure('Position',[100,100,800,600],'Color','w');grid on; box off;
for i = 1:end_steps
    plot([0; xm; L],[phi_0; phi(:,i); phi_L], ...
        '-.','LineWidth',2, ...
        'DisplayName',sprintf('t = %.4f s',dt*i));

    ax=gca;
    ax.FontName='default';ax.FontSize=12;
    ax.FontWeight="bold";ax.FontAngle="italic";
    title(sprintf('FVM solution of 1D unsteady scalar transport'),'FontAngle','normal');
    xlabel('x','FontSize',13.2);ylabel('\boldmath$\phi$','Interpreter','latex','FontSize',13.2);
    xticks(0:L/5:L);yticks(phi_0:(phi_L-phi_0)/5:phi_L);

    %legend('Location','northwest');
    %pause(0.01);
    drawnow;

end
%% animation (smart updating)
f3=figure('Position',[100,100,800,600],'Color','w');grid on; box off;
pt = plot([0; xm; L],[phi_0; phi(:,1); phi_L], ...
    '-.','LineWidth',2, ...
    'DisplayName',sprintf('t = %.4f s',dt*1));

ax=gca;
ax.FontName='default';ax.FontSize=12;
ax.FontWeight="bold";ax.FontAngle="italic";
title(sprintf('FVM solution of 1D unsteady scalar transport'),'FontAngle','normal');
xlabel('x','FontSize',13.2);ylabel('\boldmath$\phi$','Interpreter','latex','FontSize',13.2);
xticks(0:L/5:L);yticks(phi_0:(phi_L-phi_0)/5:phi_L);

for i = 1:end_steps
    set(pt,'YData',[phi_0; phi(:,i); phi_L],'DisplayName',sprintf('t = %.4f s',dt*i));
    %legend('Location','northwest');
    drawnow;
end