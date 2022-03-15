clear all
set(0, 'DefaultAxesFontWeight', 'bold', ...
    'DefaultAxesFontSize', 12, ...
    'DefaultAxesFontAngle', 'normal', ... % Not sure the difference here
    'DefaultAxesFontWeight', 'bold', ... % Not sure the difference here
    'DefaultAxesTitleFontWeight', 'bold', ...
    'DefaultAxesTitleFontSizeMultiplier', 1) ;
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesLineWidth', 2);


figure
set(gcf,'Position',[0, 0, 1100, 1000])

%%
run_chen_depth_n;


ax(3)=subplot(2,3,1);
hold on
plot(xdepth,Rnet1,'linewidth',2)
ylabel('R_p^{CaCO_3} (mol/m^2/s)')
xlabel('Depth (m)')
set(gca,'yscale','log')
box on
ax = gca;
ax.LineWidth = 1.5;
axis([1000,5000,-inf,inf])
ylim([0.5e-6 2e-6])
title(['(g)'],'Fontsize',14);
axis square

ax(3)=subplot(2,3,2);
hold on
plot(xdepth,Kp,'linewidth',2)
xlabel('Depth (m)')
ylabel('K_{Sr}')
box on
ax = gca;
ax.LineWidth = 1.5;
axis([1000,5000,-inf,inf])
ylim([0.13 0.18])
title(['(h)'],'Fontsize',14);
axis square


ax(3)=subplot(2,3,3);
hold on
plot(xdepth,xSr./xCa,'linewidth',2)
xlabel('Depth (m)')
ylabel('Sr/Ca_{cf}')
box on
ax = gca;
ax.LineWidth = 1.5;
axis([1000,5000,8.1e-3,8.6e-3])
ylim([8.5e-3 8.7e-3])
title(['(i)'],'Fontsize',14);
axis square
%%
run_chen_depth_a;


ax(3)=subplot(2,3,1);
hold on
plot(xdepth,Rnet1,'linewidth',2)

ax(3)=subplot(2,3,2);
hold on
plot(xdepth,Kp,'linewidth',2)

ax(3)=subplot(2,3,3);
hold on
plot(xdepth,xSr./xCa,'linewidth',2)
%%
run_chen_depth_i;


ax(3)=subplot(2,3,1);
hold on
plot(xdepth,Rnet1,'linewidth',2)

ax(3)=subplot(2,3,2);
hold on
plot(xdepth,Kp,'linewidth',2)

ax(3)=subplot(2,3,3);
hold on
plot(xdepth,xSr./xCa,'linewidth',2)
%%
run_chen_depth_p;
load carbonate_parameters;

ax(3)=subplot(2,3,1);
hold on
plot(xdepth,Rnet1,'linewidth',2)
legend('Norwegian Sea','Atlantic Ocean','Indian Ocean','Pacific Ocean','location','southwest','fontsize',10);

ax(3)=subplot(2,3,2);
hold on
plot(xdepth,Kp,'linewidth',2)
set(gca,'ColorOrderIndex',1)
scatter(WD_yu_n,ksr_yu_n,'o')
scatter(WD_yu_a,ksr_yu_a,'s')
scatter(WD_yu_i,ksr_yu_i,'d')
scatter(WD_yu_p,ksr_yu_p,'p')

ax(3)=subplot(2,3,3);
hold on
plot(xdepth,xSr./xCa,'linewidth',2)
xx = 1000:1000:5000;
yy = 0.00854+zeros(5,1);
plot(xx,yy,'--','color','k')
text(2000,0.00855,' seawater','FontName','Times New Roman','fontweight','bold','fontsize',12)


set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14,'Fontweight','bold','LineWidth', 2)
