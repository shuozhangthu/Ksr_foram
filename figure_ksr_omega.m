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
set(gcf,'Position',[100, 0, 700,600])

%%
run_chen_depth_n;
hold on
plot(xomega_sw,Kp,'linewidth',2)
xlabel('\Omega_{sw}','fontsize',18,'interpreter','tex')
ylabel('K_{Sr}','fontsize',18,'interpreter','tex')


box on
ax = gca;
set(gca,'yscale','log')
box on
ax = gca;
ax.LineWidth = 1.5;

%%
run_chen_depth_a;


hold on
plot(xomega_sw,Kp,'linewidth',2)

%%
run_chen_depth_i;


hold on
plot(xomega_sw,Kp,'linewidth',2)

%%
run_chen_depth_p;
load carbonate_parameters;

hold on
plot(xomega_sw,Kp,'linewidth',2)
set(gca,'ColorOrderIndex',1)
scatter(omega_yu_n,ksr_yu_n,'o','linewidth',1);
scatter(omega_yu_a,ksr_yu_a,'s','linewidth',1);
scatter(omega_yu_i,ksr_yu_i,'d','linewidth',1);
scatter(omega_yu_p,ksr_yu_p,'p','linewidth',1);
legend('Norwegian Sea','Atlantic Ocean','Indian Ocean','Pacific Ocean','location','southwest','fontsize',12);
xlim([0.6 2.2])
ylim([0.13 0.18])
set(gca,'YTick',0.13:0.01:0.18);

set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14,'Fontweight','bold','LineWidth', 2)