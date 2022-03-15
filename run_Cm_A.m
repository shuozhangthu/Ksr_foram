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
load carbonate_parameters_Cm 
%%
xdepth = [1614:100:5000]';
xpH_sw = -4e-5*xdepth+8.0991;
xDIC_sw= -0.0085*xdepth+2176.6;
xTemp = 4e-7*xdepth.*xdepth-0.0029*xdepth+8.2168;
P=1+xdepth/10;
xlnK_AB_A= log(10^(-8.3889)) +44.4236*P/83.15/276.22-13.3178*P.*P*0.001/2/83.15/276.22;
xlnK_hco3_A= log(10^(-10.5822)) +16.8052*P/83.15/276.22-4.87*P.*P*0.001/2/83.15/276.22;
xlnK_MB_A= xlnK_AB_A+ 2.1817 +((P-1)*5.1*0.01/8.315/276.22);

xK_AB_sw=exp(xlnK_AB_A);
xK_MB_sw=exp(xlnK_MB_A);
xK_hco3_sw=exp(xlnK_hco3_A);

 
 
xCO2 = zeros(length(xdepth),1);
xCO3 = zeros(length(xdepth),1);
xHCO3 = zeros(length(xdepth),1);
xomega = zeros(length(xdepth),1);
xDIC = zeros(length(xdepth),1);
xCa = zeros(length(xdepth),1);
xSr = zeros(length(xdepth),1);
xpH = zeros(length(xdepth),1);
xTA = zeros(length(xdepth),1);
xK_AB=zeros(length(xdepth),1);
xK_MB=zeros(length(xdepth),1);
xK_hco3=zeros(length(xdepth),1);


for ii = 1:length(xdepth)
    
    [o,time,out]=function_Cm_A(xdepth(ii),xpH_sw(ii),xDIC_sw(ii),xTemp(ii),xK_AB_sw(ii),xK_MB_sw(ii),xK_hco3_sw(ii));
    
    DIC = zeros(length(out(:,1)),1);
    pH = zeros(length(out(:,1)),1);
    H = zeros(length(out(:,1)),1);
    chi = zeros(length(out(:,1)),1);
    omega = zeros(length(out(:,1)),1);
    CO2 = zeros(length(out(:,1)),1);
    Ca = zeros(length(out(:,1)),1);
    Sr = zeros(length(out(:,1)),1);
    JCaCO3 = zeros(length(out(:,1)),1);
    TA = zeros(length(out(:,1)),1);
    CO3 = zeros(length(out(:,1)),1);
    HCO3 = zeros(length(out(:,1)),1);
    E2666 = zeros(length(out(:,1)),1);
    K_AB=zeros(length(out(:,1)),1);
    K_MB=zeros(length(out(:,1)),1);
    K_hco3=zeros(length(out(:,1)),1);
    
    
    
    for i = 1:length(out(:,1))
        DIC(i) = out(i,1)+out(i,2);
        a = 1;
        KN = 0;
        NT = 0;
        K1 = o.K1;
        K2 = o.K2;
        Kw = o.Kw;
        Ksp = o.Ksp;
        TA(i) = out(i,4);
        b = KN+TA(i)+K1;
        c = TA(i)*KN-KN*NT-Kw+TA(i)*K1+KN*K1+K1*K2...
            -DIC(i)*K1;
        d = TA(i)*KN*K1+TA(i)*K1*K2-Kw*KN-KN*NT*K1-Kw*K1+KN*K1*K2...
            -DIC(i)*KN*K1-2*DIC(i)*K1*K2;
        e = TA(i)*KN*K1*K2-Kw*KN*K1-KN*NT*K1*K2-Kw*K1*K2...
            -DIC(i)*KN*2*K1*K2;
        f = -K1*K2*Kw*KN;
        p = [a b c d e f];
        r = roots(p);
        H(i) = max(real(r));
        pH(i) = -log10(H(i));
        chi(i) = 1/(1+K2/H(i));
        Ca(i) = out(i,3);
        Sr(i) = out(i,5);
        CO2(i) = out(i,1);
        HCO3(i) = out(i,2)/(1+K2/H(i));
        CO3(i) = out(i,2)-HCO3(i);
        E2666(i) = out(i,2);
        omega(i) = Ca(i)*CO3(i)/Ksp;
    end
    
    xCO2(ii) = CO2(end);
    xCO3(ii) = CO3(end);
    xHCO3(ii) = HCO3(i);
    xDIC(ii) = DIC(end);
    xCa(ii) = Ca(end);
    xSr(ii) = Sr(end);
    xomega(ii) = omega(end);
    xpH(ii) = pH(end);
    xTA(ii) = TA(end);
    xomega_sw(ii)=o.omega;
    xDIC_sw(ii)=o.DIC;
    xpH_sw(ii)=o.pH;
    xCO3_sw(ii)=o.CO3;
    xK_AB_sw(ii)=o.K_AB;
    xK_MB_sw(ii)=o.K_MB;
    xK_hco3_sw(ii)=o.K_hco3;
    xCO2_sw(ii)=o.CO2;
    xHCO3_sw(ii)=o.HCO3;
    
    
A=o.A;
M=o.M;
K_AB=xK_AB_sw(ii);
K_MB=xK_MB_sw(ii);
a0=o.a0;
boltz=o.boltz;% Boltzmann constant
TK=o.TK;
a_new=o.a_new;
b_new=o.b_new;
h=o.h;
d_new=o.d_new;
alpha=o.alpha;
gamma=o.gamma;
epsilon=o.epsilon;
phi1=-log10(xK_hco3_sw(ii));
theta=10^(8.6-xpH(ii));
phi=10^(phi1-xpH(ii));
    
k_A1=o.k_A1;
k_A2=k_A1;
k_B1=2*(1+theta)*k_A1/(1+phi);
k_B2=k_B1;
k_M1=0.3*k_A1;
k_M2=k_M1;
k_BM1=k_B1;
k_BM2=k_BM1;

k_A=k_A1+k_A2*theta;
k_B=k_B1+k_B2*phi;
k_M=k_M1+k_M2*theta;
k_BM=k_BM1+k_BM2*phi;

v_A1=o.v_A1;
v_A2=v_A1;
v_A=v_A1+v_A2;

v_M1=o.v_M1;
v_M2=v_M1;
v_M=v_M1+v_M2;
    
    
    
    omega_new=xomega(ii);
    
    B=omega_new*K_AB/A;
    
    P_A=0.5;
    P_M=0.01;
    P_B=(1-P_A-P_M)/(1+theta);
    x=0;
    
    for i=1:20
        
        v_B=exp(a0*x^2)*K_AB*k_A*k_B/v_A;
        v_BM=exp(a0*(1-x)^2)*K_MB*k_M*k_BM/v_M;
        
        [P_A,P_M,P_B]=solveP(k_A,k_B,k_M,k_BM,v_A,v_B,v_M,v_BM,P_A,P_M,A,B,M,theta);

        u_A=k_A*A*P_B-v_A*P_A;
        u_M=k_M*M*P_B-v_M*P_M;
        
        rx=u_M/u_A;
        x=rx/(1+rx);
        
        Kp(ii)=rx/M*A;
        
        u_BA=k_B*B*P_A-v_B*P_B*(1-x);
        u_BM=k_BM*B*P_M-v_BM*P_B*x;
        
        u_net=u_A+u_M+u_BA+u_BM;
        
        Omega_A=2+(v_A*exp(2*epsilon/boltz/TK)/k_B/B+v_A*exp(2*epsilon/boltz/TK)*v_B/k_A/A/k_B/B)/(1-v_A*v_B/k_A/A/k_B/B);
        Omega_B=2+(v_B*exp(2*epsilon/boltz/TK)/k_A/A+v_A*exp(2*epsilon/boltz/TK)*v_B/k_A/A/k_B/B)/(1-v_A*v_B/k_A/A/k_B/B);
        
        C=2/u_net*(k_A*A/Omega_A+k_B*B/Omega_B);
        
        rho=((C^2+4*C)^(0.5)-C)/2;
        
        v_st=rho*u_net*b_new;
        sigma=abs(log(omega_new));
        beta_st=v_st/sigma;
        Delta_G=1/8*h*a_new*b_new*gamma*gamma/h/boltz/TK/sigma;
        I=beta_st*h/h/a_new/b_new*(4*h*sigma/3.14/h/a_new/b_new)^0.5*exp(-Delta_G/boltz/TK);
        
        Rnet1(ii) = 1.137*h*(I*v_st^2)^(1/3)*d_new; % Calcite growth rate
        
        y0=4*2*h*h*a_new*alpha/(1.3807e-23)/TK/sigma;
        Rnet(ii)=rho*u_net*h*b_new*d_new/2/y0;
            

    end
    
    
end
% %
set(gcf,'DefaultAxesLineWidth',1.2)
set(gcf,'DefaultLineLineWidth',1.2)
%%
subplot(2,3,1)
hold on
xlabel('Depth (m)');
ylabel('pH_{sw}');
box on
ax = gca;
ax.LineWidth = 1.5;
plot(xdepth,xpH_sw,'k','linewidth',2);
axis([1000,5000,-inf,inf])
title(['(a)'],'Fontsize',14);
ylim([7.9 8.05])


subplot(2,3,2)
hold on
xlabel('Depth (m)');
ylabel('DIC_{sw} (\mumol/kg)');
box on
ax = gca;
ax.LineWidth = 1.5;
plot(xdepth,xDIC_sw*1000000,'k','linewidth',2);
axis([1000,5000,-inf,inf])
title(['(b)'],'Fontsize',14);
ylim([2100 2400])


subplot(2,3,3)
hold on
xlabel('Depth (m)');
ylabel('\Omega_{sw}');
box on
ax = gca;
ax.LineWidth = 1.5;
plot(xdepth,xomega_sw,'k','linewidth',2)
axis([1000,5000,-inf,inf])
ylim([1 2.5])
title(['(c)'],'Fontsize',14);


subplot(2,3,4);
hold on
plot(xdepth,xpH,'k','linewidth',2)
ylabel('pH_{cf}')
xlabel('Depth (m)');
box on
ax = gca;
ax.LineWidth = 1.5;
axis([1000,5000,-inf,inf])
ylim([9.12 9.18])
title(['(d)'],'Fontsize',14);

ax(3)=subplot(2,3,5);
hold on
plot(xdepth,xDIC*1000000,'k','linewidth',2)
ylabel('DIC_{cf} (\mumol/kg)')
xlabel('Depth (m)');
box on
ax = gca;
ax.LineWidth = 1.5;
axis([1000,5000,-inf,inf])
ylim([2020 2100])
title(['(e)'],'Fontsize',14);

ax(3)=subplot(2,3,6);
hold on
plot(xdepth,xomega,'k','linewidth',2)
ylabel('\Omega_{cf}')
xlabel('Depth (m)');
box on
ax = gca;
ax.LineWidth = 1.5;
axis([1000,5000,-inf,inf])
ylim([8 16])
title(['(f)'],'Fontsize',14);

load carbonate_parameters_Cm 
subplot(2,3,1)
hold on
set(gca,'ColorOrderIndex',1)
scatter(WD_A,pH_A,'o','k')


subplot(2,3,2)
hold on
set(gca,'ColorOrderIndex',1)
scatter(WD_A,DIC_A,'o','k');


subplot(2,3,3)
hold on
set(gca,'ColorOrderIndex',1)
scatter(WD_A,omega_cm_A,'o','k')
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14,'Fontweight','bold','LineWidth', 2)


figure
set(gcf,'Position',[0, 0, 1100, 1000])

load carbonate_parameters_Cm 

ax(3)=subplot(2,3,1);
hold on
plot(xdepth,Rnet1,'k','linewidth',2)
ylabel('R_p^{CaCO_3} (mol/m^2/s)')
xlabel('Water depth (m)')
set(gca,'yscale','log')
box on
ax = gca;
ax.LineWidth = 1.5;
title(['(g)'],'Fontsize',14);
axis([1000,5000,-inf,inf])
ylim([3e-7 15e-7])

ax(3)=subplot(2,3,2);
hold on
plot(xdepth,Kp,'k','linewidth',2)
xlabel('Water depth (m)')
ylabel('K_{Sr}')
box on
ax = gca;
ax.LineWidth = 1.5;
title(['(h)'],'Fontsize',14);
axis([1000,5000,-inf,inf])
ylim([0.12 0.16])

ax(3)=subplot(2,3,3);
hold on
plot(xomega_sw,Kp,'k','linewidth',2)
xlabel('\Omega_{sw}')
ylabel('K_{Sr}')
box on
ax = gca;
ax.LineWidth = 1.5;
title(['(i)'],'Fontsize',14);
ylim([0.12 0.16])
xlim([1 2.5])


subplot(2,3,2)
hold on
set(gca,'ColorOrderIndex',1)
scatter(WD_A,DSr_cm_A,'o','k')

subplot(2,3,3)
hold on
set(gca,'ColorOrderIndex',1)
scatter(omega_cm_A,DSr_cm_A,'o','k');


set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14,'Fontweight','bold','LineWidth', 2)