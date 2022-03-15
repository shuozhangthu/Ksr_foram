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

subplot(2,3,1)
hold on
xlabel('Depth (m)');
ylabel('pH_{sw}');
box on
ax = gca;
ax.LineWidth = 1.5;
plot(xdepth,xpH_sw);
axis([1000,5000,-inf,inf])
ylim([7.7 8.2])
title(['(a)'],'Fontsize',14);
axis square

subplot(2,3,2)
hold on
xlabel('Depth (m)');
ylabel('DIC_{sw} (\mumol/kg)');
box on
ax = gca;
ax.LineWidth = 1.5;
plot(xdepth,xDIC_sw*1000000);
axis([1000,5000,-inf,inf])
ylim([2050 2350])
title(['(b)'],'Fontsize',14);
axis square

subplot(2,3,3)
hold on
xlabel('Depth (m)');
ylabel('\Omega_{sw}');
box on
ax = gca;
ax.LineWidth = 1.5;
plot(xdepth,xomega_sw,'linewidth',2)
axis([1000,5000,-inf,inf])
ylim([0.5 2.5])
title(['(c)'],'Fontsize',14);
axis square

subplot(2,3,4);
hold on
plot(xdepth,xpH,'linewidth',2)
ylabel('pH_{cf}')
xlabel('Depth (m)');
box on
ax = gca;
ax.LineWidth = 1.5;
axis([1000,5000,-inf,inf])
ylim([9.2 9.5])
title(['(d)'],'Fontsize',14);
axis square

ax(3)=subplot(2,3,5);
hold on
plot(xdepth,xDIC*1000000,'linewidth',2)
ylabel('DIC_{cf} (\mumol/kg)')
xlabel('Depth (m)');
box on
ax = gca;
ax.LineWidth = 1.5;
axis([1000,5000,-inf,inf])
ylim([1959 2250])
title(['(e)'],'Fontsize',14);
axis square

ax(3)=subplot(2,3,6);
hold on
plot(xdepth,xomega,'linewidth',2)
ylabel('\Omega_{cf}')
xlabel('Depth (m)');
box on
ax = gca;
ax.LineWidth = 1.5;
axis([1000,5000,-inf,inf])
ylim([10 20])
title(['(f)'],'Fontsize',14);
axis square

%%
run_chen_depth_a;

subplot(2,3,1)
hold on
plot(xdepth,xpH_sw);

subplot(2,3,2)
hold on
plot(xdepth,xDIC_sw*1000000);

subplot(2,3,3)
plot(xdepth,xomega_sw,'linewidth',2)

subplot(2,3,4);
hold on
plot(xdepth,xpH,'linewidth',2)

ax(3)=subplot(2,3,5);
hold on
plot(xdepth,xDIC*1000000,'linewidth',2)

ax(3)=subplot(2,3,6);
hold on
plot(xdepth,xomega,'linewidth',2)

%%
run_chen_depth_i;

subplot(2,3,1)
hold on
plot(xdepth,xpH_sw);

subplot(2,3,2)
hold on
plot(xdepth,xDIC_sw*1000000);

subplot(2,3,3)
plot(xdepth,xomega_sw,'linewidth',2)

subplot(2,3,4);
hold on
plot(xdepth,xpH,'linewidth',2)

ax(3)=subplot(2,3,5);
hold on
plot(xdepth,xDIC*1000000,'linewidth',2)

ax(3)=subplot(2,3,6);
hold on
plot(xdepth,xomega,'linewidth',2)

%%
run_chen_depth_p;
load carbonate_parameters;

subplot(2,3,1)
hold on
plot(xdepth,xpH_sw);
set(gca,'ColorOrderIndex',1)
scatter(WD_yu_n,pH_yu_n,'o')
scatter(WD_yu_a,pH_yu_a,'s')
scatter(WD_yu_i,pH_yu_i,'d')
scatter(WD_yu_p,pH_yu_p,'p')

subplot(2,3,2)
hold on
plot(xdepth,xDIC_sw*1000000);
set(gca,'ColorOrderIndex',1)
scatter(WD_yu_n,DIC_yu_n,'o');
scatter(WD_yu_a,DIC_yu_a,'s');
scatter(WD_yu_i,DIC_yu_i,'d');
scatter(WD_yu_p,DIC_yu_p,'p');

subplot(2,3,3)
hold on
plot(xdepth,xomega_sw,'linewidth',2)

set(gca,'ColorOrderIndex',1)
scatter(WD_yu_n,omega_yu_n,'o')
scatter(WD_yu_a,omega_yu_a,'s')
scatter(WD_yu_i,omega_yu_i,'d')
scatter(WD_yu_p,omega_yu_p,'p')


subplot(2,3,4);
hold on
plot(xdepth,xpH,'linewidth',2)

ax(3)=subplot(2,3,5);
hold on
plot(xdepth,xDIC*1000000,'linewidth',2)

ax(3)=subplot(2,3,6);
hold on
plot(xdepth,xomega,'linewidth',2)

set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14,'Fontweight','bold','LineWidth', 2)


