clear
%kA
kA_xdata = ...
 [1/278.15	1/298.15	1/313.15];
kA_ydata = ...
 [10^(5.525)	10^(5.5796)	10^(5.6286)];

%va
va_xdata = ...
 [1/278.15	1/298.15	1/313.15];
va_ydata = ...
[10^(0.308)	10^(0.3539)	10^(0.4611)];

%vm
vm_xdata = ...
[1/278.15	1/298.15	1/313.15];
vm_ydata = ...
 [10^(0.5953)	10^(0.8342)	10^(1.1532)];
figure 
hold on 
ka_x=[3.1e-3:0.00001:3.6e-3];
ka_y=2.824e6*exp(-594.6*ka_x);
scatter(kA_xdata,kA_ydata)
scatter(va_xdata,va_ydata)
scatter(vm_xdata,vm_ydata)

set(gca,'ColorOrderIndex',1)

plot(ka_x,ka_y,':','linewidth',2)

va_x=[3.1e-3:0.00001:3.6e-3];
va_y=50.55*exp(-905.5*va_x);
plot(va_x,va_y,'-','linewidth',2)

vm_x=[3.1e-3:0.00001:3.6e-3];
vm_y=2.085e6*exp(-3731*vm_x);
plot(vm_x,vm_y,'--','linewidth',2)

xlabel('1/T (K^{-1})')
ylabel('ν (s^{-1}) or k (s^{-1}M^{-1}) ')
set(gca,'yscale','log')
box on
ax = gca;
ax.LineWidth = 1.5;
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14,'Fontweight','bold','LineWidth', 2)
xlim([3.1e-3,3.6e-3]);

legend('k_{A1}','ν_{A1}','ν_{M1}','location','best')

T5=278.15;
T25=298.15;
T40=313.15;
T_Norwegian=272.26;
T_Atlantic=276.01;
T_Indian=274.51;
T_Pacific=275.05;
T_A=3.07+273.15;
T_I=1.33+273.15;
T_P=1.54+273.15;
%Kcaco3
logKcaco3_T5=-171.9065-0.077993*T5+2839.3/T5+71.595*log10(T5);
logKcaco3_T25=-171.9065-0.077993*T25+2839.3/T25+71.595*log10(T25);
logKcaco3_T40=-171.9065-0.077993*T40+2839.3/T40+71.595*log10(T40);
logKcaco3_T_Norwegian=-171.9065-0.077993*T_Norwegian+2839.3/T_Norwegian+71.595*log10(T_Norwegian);
logKcaco3_T_Atlantic=-171.9065-0.077993*T_Atlantic+2839.3/T_Atlantic+71.595*log10(T_Atlantic);
logKcaco3_T_Indian=-171.9065-0.077993*T_Indian+2839.3/T_Indian+71.595*log10(T_Indian);
logKcaco3_T_Pacific=-171.9065-0.077993*T_Pacific+2839.3/T_Pacific+71.595*log10(T_Pacific);

logKcaco3_T_A=-171.9065-0.077993*T_A+2839.3/T_A+71.595*log10(T_A);
logKcaco3_T_I=-171.9065-0.077993*T_I+2839.3/T_I+71.595*log10(T_I);
logKcaco3_T_P=-171.9065-0.077993*T_P+2839.3/T_P+71.595*log10(T_P);



Kcaco3_T5=10^logKcaco3_T5;
Kcaco3_T25=10^logKcaco3_T25;
Kcaco3_T40=10^logKcaco3_T40;
Kcaco3_T_Norwegian=10^logKcaco3_T_Norwegian;
Kcaco3_T_Atlantic=10^logKcaco3_T_Atlantic;
Kcaco3_T_Indian=10^logKcaco3_T_Indian;
Kcaco3_T_Pacific=10^logKcaco3_T_Pacific;

Kcaco3_T_A=10^logKcaco3_T_A;
Kcaco3_T_I=10^logKcaco3_T_I;
Kcaco3_T_P=10^logKcaco3_T_P;


%Khco3
logKhco3_T5=-107.8871-0.03252894*T5+5151.79/T5+38.92561*log10(T5)-563713.9/T5/T5;
logKhco3_T25=-107.8871-0.03252894*T25+5151.79/T25+38.92561*log10(T25)-563713.9/T25/T25;
logKhco3_T40=-107.8871-0.03252894*T40+5151.79/T40+38.92561*log10(T40)-563713.9/T40/T40;
logKhco3_T_Norwegian=-107.8871-0.03252894*T_Norwegian+5151.79/T_Norwegian+38.92561*log10(T_Norwegian)-563713.9/T_Norwegian/T_Norwegian;
logKhco3_T_Atlantic=-107.8871-0.03252894*T_Atlantic+5151.79/T_Atlantic+38.92561*log10(T_Atlantic)-563713.9/T_Atlantic/T_Atlantic;
logKhco3_T_Indian=-107.8871-0.03252894*T_Indian+5151.79/T_Indian+38.92561*log10(T_Indian)-563713.9/T_Indian/T_Indian;
logKhco3_T_Pacific=-107.8871-0.03252894*T_Pacific+5151.79/T_Pacific+38.92561*log10(T_Pacific)-563713.9/T_Pacific/T_Pacific;

logKhco3_T_A=-107.8871-0.03252894*T_A+5151.79/T_A+38.92561*log10(T_A)-563713.9/T_A/T_A;
logKhco3_T_I=-107.8871-0.03252894*T_I+5151.79/T_I+38.92561*log10(T_I)-563713.9/T_I/T_I;
logKhco3_T_P=-107.8871-0.03252894*T_P+5151.79/T_P+38.92561*log10(T_P)-563713.9/T_P/T_P;


Khco3_T_Norwegian=10^(logKhco3_T_Norwegian);
Khco3_T_Atlantic=10^(logKhco3_T_Atlantic);
Khco3_T_Indian=10^(logKhco3_T_Indian);
Khco3_T_Pacific=10^(logKhco3_T_Pacific);

Khco3_T_A=10^(logKhco3_T_A);
Khco3_T_I=10^(logKhco3_T_I);
Khco3_T_P=10^(logKhco3_T_P);

ksr_a=log( 3.3113e-09)-log(2.5*10^(-8));
ksr_a=exp(ksr_a);

%Ksrco3
ksr_T5=ksr_a*exp(5000/8.315*(1/298.15-1/T5)); 
ksr_T25=ksr_a*exp(5000/8.315*(1/298.15-1/T25));
ksr_T40=ksr_a*exp(5000/8.315*(1/298.15-1/T40));
ksr_T_Norwegian=ksr_a*exp(5000/8.315*(1/298.15-1/T_Norwegian)); 
ksr_T_Atlantic=ksr_a*exp(5000/8.315*(1/298.15-1/T_Atlantic)); 
ksr_T_Indian=ksr_a*exp(5000/8.315*(1/298.15-1/T_Indian)); 
ksr_T_Pacific=ksr_a*exp(5000/8.315*(1/298.15-1/T_Pacific)); 

ksr_T_A=ksr_a*exp(5000/8.315*(1/298.15-1/T_A)); 
ksr_T_I=ksr_a*exp(5000/8.315*(1/298.15-1/T_I)); 
ksr_T_P=ksr_a*exp(5000/8.315*(1/298.15-1/T_P)); 


lnKsrco3_T5=log(Kcaco3_T5)-log(ksr_T5);
lnKsrco3_T25=log(Kcaco3_T25)-log(ksr_T25);
lnKsrco3_T40=log(Kcaco3_T40)-log(ksr_T40);
lnKsrco3_T_Norwegian=log(Kcaco3_T_Norwegian)-log(ksr_T_Norwegian);
lnKsrco3_T_Atlantic=log(Kcaco3_T_Atlantic)-log(ksr_T_Atlantic);
lnKsrco3_T_Indian=log(Kcaco3_T_Indian)-log(ksr_T_Indian);
lnKsrco3_T_Pacific=log(Kcaco3_T_Pacific)-log(ksr_T_Pacific);

lnKsrco3_T_A=log(Kcaco3_T_A)-log(ksr_T_A);
lnKsrco3_T_I=log(Kcaco3_T_A)-log(ksr_T_I);
lnKsrco3_T_P=log(Kcaco3_T_A)-log(ksr_T_P);




Ksrco3_T5=exp(lnKsrco3_T5);
Ksrco3_T25=exp(lnKsrco3_T25);
Ksrco3_T40=exp(lnKsrco3_T40);
Ksrco3_T_Norwegian=exp(lnKsrco3_T_Norwegian);
Ksrco3_T_Atlantic=exp(lnKsrco3_T_Atlantic);
Ksrco3_T_Indian=exp(lnKsrco3_T_Indian);
Ksrco3_T_Pacific=exp(lnKsrco3_T_Pacific);

Ksrco3_T_A=exp(lnKsrco3_T_A);
Ksrco3_T_I=exp(lnKsrco3_T_I);
Ksrco3_T_P=exp(lnKsrco3_T_P);


t_Norwegian=-0.89;
t_Atlantic=2.86;
t_Indian=1.36;
t_Pacific=1.9;
t_A=3.07;
t_I=1.33;
t_P=1.54;


Vca_t_Norwegian=-15.93-0.0196*t_Norwegian+0.668*0.001*t_Norwegian*t_Norwegian; 
Vca_t_Atlantic=-15.93-0.0196*t_Atlantic+0.668*0.001*t_Atlantic*t_Atlantic;
Vca_t_Indian=-15.93-0.0196*t_Indian+0.668*0.001*t_Indian*t_Indian;
Vca_t_Pacific=-15.93-0.0196*t_Pacific+0.668*0.001*t_Pacific*t_Pacific;
Vca_t_A=-15.93-0.0196*t_A+0.668*0.001*t_A*t_A;
Vca_t_I=-15.93-0.0196*t_I+0.668*0.001*t_I*t_I;
Vca_t_P=-15.93-0.0196*t_P+0.668*0.001*t_P*t_P;


             
Vco3_t_Norwegian=7.4+0.3725*t_Norwegian-5.653*0.001*t_Norwegian*t_Norwegian;
Vco3_t_Atlantic=7.4+0.3725*t_Atlantic-5.653*0.001*t_Atlantic*t_Atlantic;
Vco3_t_Indian=7.4+0.3725*t_Indian-5.653*0.001*t_Indian*t_Indian;
Vco3_t_Pacific=7.4+0.3725*t_Pacific-5.653*0.001*t_Pacific*t_Pacific;
Vco3_t_A=7.4+0.3725*t_A-5.653*0.001*t_A*t_A;
Vco3_t_I=7.4+0.3725*t_I-5.653*0.001*t_I*t_I;
Vco3_t_P=7.4+0.3725*t_P-5.653*0.001*t_P*t_P;



delta_V_Norwegian=Vca_t_Norwegian+Vco3_t_Norwegian-36.93;
delta_V_Atlantic=Vca_t_Atlantic+Vco3_t_Atlantic-36.93;
delta_V_Indian=Vca_t_Indian+Vco3_t_Indian-36.93;
delta_V_Pacific=Vca_t_Pacific+Vco3_t_Pacific-36.93;
delta_V_A=Vca_t_A+Vco3_t_A-36.93;
delta_V_I=Vca_t_I+Vco3_t_I-36.93;
delta_V_P=Vca_t_P+Vco3_t_P-36.93;


%K
Kca_t_Norwegian=-5.82-0.0428*t_Norwegian; 
Kca_t_Atlantic=-5.82-0.0428*t_Atlantic;
Kca_t_Indian=-5.82-0.0428*t_Indian;
Kca_t_Pacific=-5.82-0.0428*t_Pacific;
Kca_t_A=-5.82-0.0428*t_A;
Kca_t_I=-5.82-0.0428*t_I;
Kca_t_P=-5.82-0.0428*t_P;



Kco3_t_Norwegian=-7.88+0.1673*t_Norwegian;
Kco3_t_Atlantic=-7.88+0.1673*t_Atlantic;
Kco3_t_Indian=-7.88+0.1673*t_Indian;
Kco3_t_Pacific=-7.88+0.1673*t_Pacific;
Kco3_t_A=-7.88+0.1673*t_A;
Kco3_t_I=-7.88+0.1673*t_I;
Kco3_t_P=-7.88+0.1673*t_P;



delta_K_Norwegian=Kca_t_Norwegian+Kco3_t_Norwegian-0;
delta_K_Atlantic=Kca_t_Atlantic+Kco3_t_Atlantic-0;
delta_K_Indian=Kca_t_Indian+Kco3_t_Indian-0;
delta_K_Pacific=Kca_t_Pacific+Kco3_t_Pacific-0;
delta_K_A=Kca_t_A+Kco3_t_A-0;
delta_K_I=Kca_t_I+Kco3_t_I-0;
delta_K_P=Kca_t_P+Kco3_t_P-0;
%△V

Vhco3_t_Norwegian=24.72+0.1967*t_Norwegian-3.006*0.001*t_Norwegian*t_Norwegian; 
Vhco3_t_Atlantic=24.72+0.1967*t_Atlantic-3.006*0.001*t_Atlantic*t_Atlantic;
Vhco3_t_Indian=24.72+0.1967*t_Indian-3.006*0.001*t_Indian*t_Indian;
Vhco3_t_Pacific=24.72+0.1967*t_Pacific-3.006*0.001*t_Pacific*t_Pacific;
Vhco3_t_A=24.72+0.1967*t_A-3.006*0.001*t_A*t_A;
Vhco3_t_I=24.72+0.1967*t_I-3.006*0.001*t_I*t_I;
Vhco3_t_P=24.72+0.1967*t_P-3.006*0.001*t_P*t_P;


Vco3_t_Norwegian=7.4+0.3725*t_Norwegian-5.653*0.001*t_Norwegian*t_Norwegian;
Vco3_t_Atlantic=7.4+0.3725*t_Atlantic-5.653*0.001*t_Atlantic*t_Atlantic;
Vco3_t_Indian=7.4+0.3725*t_Indian-5.653*0.001*t_Indian*t_Indian;
Vco3_t_Pacific=7.4+0.3725*t_Pacific-5.653*0.001*t_Pacific*t_Pacific;
Vco3_t_A=7.4+0.3725*t_A-5.653*0.001*t_A*t_A;
Vco3_t_I=7.4+0.3725*t_I-5.653*0.001*t_I*t_I;
Vco3_t_P=7.4+0.3725*t_P-5.653*0.001*t_P*t_P;



delta_V2_Norwegian=Vco3_t_Norwegian-Vhco3_t_Norwegian;
delta_V2_Atlantic=Vco3_t_Atlantic-Vhco3_t_Atlantic;
delta_V2_Indian=Vco3_t_Indian-Vhco3_t_Indian;
delta_V2_Pacific=Vco3_t_Pacific-Vhco3_t_Pacific;
delta_V2_A=Vco3_t_A-Vhco3_t_A;
delta_V2_I=Vco3_t_I-Vhco3_t_I;
delta_V2_P=Vco3_t_P-Vhco3_t_P;


Khco3_t_Norwegian=-3.01+0.1673*t_Norwegian; 
Khco3_t_Atlantic=-3.01+0.1673*t_Atlantic;
Khco3_t_Indian=-3.01+0.1673*t_Indian;
Khco3_t_Pacific=-3.01+0.1673*t_Pacific;
Khco3_t_A=-3.01+0.1673*t_A;
Khco3_t_I=-3.01+0.1673*t_I;
Khco3_t_P=-3.01+0.1673*t_P;



Kco3_t_Norwegian=-7.88+0.1673*t_Norwegian;
Kco3_t_Atlantic=-7.88+0.1673*t_Atlantic;
Kco3_t_Indian=-7.88+0.1673*t_Indian;
Kco3_t_Pacific=-7.88+0.1673*t_Pacific;
Kco3_t_A=-7.88+0.1673*t_A;
Kco3_t_I=-7.88+0.1673*t_I;
Kco3_t_P=-7.88+0.1673*t_P;



delta_K2_Norwegian=Kco3_t_Norwegian-Khco3_t_Norwegian;
delta_K2_Atlantic=Kco3_t_Atlantic-Khco3_t_Atlantic;
delta_K2_Indian=Kco3_t_Indian-Khco3_t_Indian;
delta_K2_Pacific=Kco3_t_Pacific-Khco3_t_Pacific;
delta_K2_A=Kco3_t_A-Khco3_t_A;
delta_K2_I=Kco3_t_I-Khco3_t_I;
delta_K2_P=Kco3_t_P-Khco3_t_P;


%%
Va1_AA=50.55*exp(-905.5/T_A);
Vm1_AA=2.085e6*exp(-3731/T_A);
Ka1_AA=2.824e6*exp(-594.6/T_A);

Va1_II=50.55*exp(-905.5/T_I);
Vm1_II=2.085e6*exp(-3731/T_I);
Ka1_II=2.824e6*exp(-594.6/T_I);

Va1_PP=50.55*exp(-905.5/T_P);
Vm1_PP=2.085e6*exp(-3731/T_P);
Ka1_PP=2.824e6*exp(-594.6/T_P);

Va1_NN=50.55*exp(-905.5/T_Norwegian);
Vm1_NN=2.085e6*exp(-3731/T_Norwegian);
Ka1_NN=2.824e6*exp(-594.6/T_Norwegian);

T_cma=3.07+273.15;
Va1_cm_A=50.55*exp(-905.5/T_cma);
Vm1_cm_A=2.085e6*exp(-3731/T_cma);
Ka1_cm_A=2.824e6*exp(-594.6/T_cma);
%log
Va1_AA=log10(Va1_AA);
Vm1_AA=log10(Vm1_AA);
Ka1_AA=log10(Ka1_AA);

Va1_II=log10(Va1_II);
Vm1_II=log10(Vm1_II);
Ka1_II=log10(Ka1_II);

Va1_PP=log10(Va1_PP);
Vm1_PP=log10(Vm1_PP);
Ka1_PP=log10(Ka1_PP);

Va1_cm_A=log10(Va1_cm_A);
Vm1_cm_A=log10(Vm1_cm_A);
Ka1_cm_A=log10(Ka1_cm_A);

Va1_NN=log10(Va1_NN);
Vm1_NN=log10(Vm1_NN);
Ka1_NN=log10(Ka1_NN);
%%
deltaG_originT_P=log(Ksrco3_T_P)-log(Kcaco3_T_P);
deltaG_originT_I=log(Ksrco3_T_I)-log(Kcaco3_T_I);
deltaG_originT_A=log(Ksrco3_T_A)-log(Kcaco3_T_A);
deltaG_originT_N=log(Ksrco3_T_Norwegian)-log(Kcaco3_T_Norwegian);