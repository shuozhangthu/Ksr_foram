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

subplot(2,3,1);
hold on
plot(xdepth,xCO2_sw,'linewidth',2)
xlabel('Depth (m)')
ylabel('[CO_{2}]_{sw}')
box on
ax = gca;
ax.LineWidth = 1.5;
axis([1000,5000,1.5e-5,3.5e-5])
title(['(a)'],'Fontsize',14);
axis square

subplot(2,3,2);
hold on
plot(xdepth,xCO3_sw,'linewidth',2)
xlabel('Depth (m)')
ylabel('[CO_{3}^{2-}]_{sw}')
box on
ax = gca;
ax.LineWidth = 1.5;
axis([1000,5000,7e-5,12e-5])
title(['(b)'],'Fontsize',14);
axis square

subplot(2,3,3);
hold on
plot(xdepth,xHCO3_sw,'linewidth',2)
xlabel('Depth (m)')
ylabel('[HCO_{3}^{-}]_{sw}')
box on
ax = gca;
ax.LineWidth = 1.5;
axis([1000,5000,-inf,inf])
ylim([1.95e-3 2.25e-3])
set(gca,'YTick',1.95e-3:0.05e-3:2.25e-3);

title(['(c)'],'Fontsize',14);
axis square

subplot(2,3,4);
hold on
plot(xdepth,xCO2,'linewidth',2)
xlabel('Depth (m)')
ylabel('[CO_{2}]_{cf}')
box on
ax = gca;
ax.LineWidth = 1.5;
axis([1000,5000,1e-6,2.5e-6])
title(['(d)'],'Fontsize',14);
axis square

subplot(2,3,5);
hold on
plot(xdepth,xCO3,'linewidth',2)
xlabel('Depth (m)')
ylabel('[CO_{3}^{2-}]_{cf}')
box on
ax = gca;
ax.LineWidth = 1.5;
axis([1000,5000,1e-3,1.2e-3])
title(['(e)'],'Fontsize',14);
axis square

subplot(2,3,6);
hold on
plot(xdepth,xHCO3,'linewidth',2)
xlabel('Depth (m)')
ylabel('[HCO_{3}^{-}_{cf}')
box on
ax = gca;
ax.LineWidth = 1.5;
axis([1000,5000,0.8e-3,1.2e-3])
title(['(f)'],'Fontsize',14);
axis square
%%
run_chen_depth_a;


subplot(2,3,1);
hold on
plot(xdepth,xCO2_sw,'linewidth',2)


subplot(2,3,2);
hold on
plot(xdepth,xCO3_sw,'linewidth',2)


subplot(2,3,3);
hold on
plot(xdepth,xHCO3_sw,'linewidth',2)



subplot(2,3,4);
hold on
plot(xdepth,xCO2,'linewidth',2)



subplot(2,3,5);
hold on
plot(xdepth,xCO3,'linewidth',2)



subplot(2,3,6);
hold on
plot(xdepth,xHCO3,'linewidth',2)

%%
run_chen_depth_i;

subplot(2,3,1);
hold on
plot(xdepth,xCO2_sw,'linewidth',2)


subplot(2,3,2);
hold on
plot(xdepth,xCO3_sw,'linewidth',2)


subplot(2,3,3);
hold on
plot(xdepth,xHCO3_sw,'linewidth',2)



subplot(2,3,4);
hold on
plot(xdepth,xCO2,'linewidth',2)



subplot(2,3,5);
hold on
plot(xdepth,xCO3,'linewidth',2)



subplot(2,3,6);
hold on
plot(xdepth,xHCO3,'linewidth',2)

%%
run_chen_depth_p;
load carbonate_parameters;

subplot(2,3,1);
hold on
plot(xdepth,xCO2_sw,'linewidth',2)


subplot(2,3,2);
hold on
plot(xdepth,xCO3_sw,'linewidth',2)


subplot(2,3,3);
hold on
plot(xdepth,xHCO3_sw,'linewidth',2)



subplot(2,3,4);
hold on
plot(xdepth,xCO2,'linewidth',2)



subplot(2,3,5);
hold on
plot(xdepth,xCO3,'linewidth',2)



subplot(2,3,6);
hold on
plot(xdepth,xHCO3,'linewidth',2)

legend('Norwegian Sea','Atlantic Ocean','Indian Ocean','Pacific Ocean','location','southwest','fontsize',10);
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14,'Fontweight','bold','LineWidth', 2)


