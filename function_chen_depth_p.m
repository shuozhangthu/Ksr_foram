function [o, time, out]=function_chen_depth_p1(depth,pH,DIC,Temp,K_AB,K_MB,K_hco3)
o=setopts(depth,pH,DIC,Temp,K_AB,K_MB,K_hco3);
[time, out]=driver(o);
end


function [time, out]=driver(o)
t=[1e-8 1e3];           % time (s)
CO20 = o.CO2;         %1 = 12C16O16O
EIC0 = o.EIC;       %2 = E12C16O16O16O
Ca0 = o.Ca;             %3 = Ca2+
Alk0 = o.TA;            %4 = Alk
Sr0=o.Sr;               %5 = Sr2+

y0=[CO20,EIC0,Ca0,Alk0,Sr0]; %IC vector
options = odeset('MaxOrder', 4, 'MaxStep', 1e3, 'RelTol', 10^-7, 'AbsTol',10^-7);
[time, out]=ode15s(@(t,y0) carbonate(t, y0, o.TK,o.kf1,o.kf4,o.kb1,o.kb4,o.K1,o.K2,o.Kw,o.Ksp,o.CO2cell,o.CO2sw,o.EICsw,o.Casw,o.Srsw,o.Alksw,o.z,o.Dcell,o.tau,o.fCa,o.FAlk,o.k_A1,o.v_A1,o.v_M1,o.K_AB,o.K_MB,o.K_hco3), t, y0, options);
end




function o=setopts(depth,pH,DIC,Temp,K_AB,K_MB,K_hco3)
%Set parameters
o.FAlk = 1.15e-5;%1.2e-6;%1e-6;              %moles/m2/s
o.fCa = 0;
o.tau = 1;%230;                %s
o.Dcell = 0;        %m/s
o.z = 10e-6;                %m
%o.Sp = 1e-2;                %m2/kgsoln (1e-2 is equivalent to z = 10 microns in Chen model)

o.TC = Temp;
o.TK = o.TC+273.15;
o.S = 34.7;
o.pH = pH;
o.DIC = DIC/1000000;

o.K_AB=K_AB;
o.K_MB=K_MB;
o.K_hco3=K_hco3;

o.Ca = 10.3e-3;%10.3e-3;
o.Sr=0.00008795;
o.krate = 10^(-0.106);
% o.krate = 0;
o.prefactor = 1000;
%ion by ion model
o.A=0.01*0.17;
o.M=8.54/1000*o.A;
o.a0=1.7;

o.boltz=1.38065E-23;% Boltzmann constant
o.TK=o.TK;
o.a_new=6.4e-10;
o.b_new=o.a_new/2;
o.h=3.1e-10;
o.d_new=27100;
o.alpha=1.41;
o.gamma=1.49e-10;
o.epsilon=10^(-19.9875);

o.phi1=-log10(o.K_hco3);
o.theta=10^(8.6-o.pH);
o.phi=10^(o.phi1-o.pH);

o.k_A1=10^( 5.5108);
o.k_A2=o.k_A1;
o.k_B1=2*(1+o.theta)*o.k_A1/(1+o.phi);
o.k_B2=o.k_B1;
o.k_M1=0.3*o.k_A1;
o.k_M2=o.k_M1;
o.k_BM1=o.k_B1;
o.k_BM2=o.k_BM1;

o.k_A=o.k_A1+o.k_A2*o.theta;
o.k_B=o.k_B1+o.k_B2*o.phi;
o.k_M=o.k_M1+o.k_M2*o.theta;
o.k_BM=o.k_BM1+o.k_BM2*o.phi;

o.v_A1=10^(0.2721);
o.v_A2=o.v_A1;
o.v_A=o.v_A1+o.v_A2;

o.v_M1=10^( 0.4203);
o.v_M2=o.v_M1;
o.v_M=o.v_M1+o.v_M2;
% 
%Equilibrium constants
P =  depth/10+1;  % bar
S=o.S;
TC=o.TC;
phflag = 0;     % 0: Total scale
% 1: Free scale
%----- choose K1 and K2: Roy or Mehrbach.
k1k2flag = 1;	% 0: Roy et al. (1993)
% 1: Mehrbach et al (1973) as
%    refit by Lueker et al. (2000) on total scale.
%----- compare output to CO2SYS? use phflag=0, k1k2flag=1
ocdflag = 1;    % 1: should be close to CO2SYS V1.1 Matlab 09/2011
% 0: else
equic;
o.K1=K1;
o.K2=K2;
o.Kw=Kw;

S=o.S;
TC=o.TC;

tmp1 = -171.9065-0.077993.*o.TK+2839.319./o.TK+71.595.*log10(o.TK);
tmp2 = +(-0.77712+0.0028426.*o.TK+178.34./o.TK).*sqrt(o.S);
tmp3 = -0.07711.*o.S+0.0041249.*o.S.^1.5;
log10Kspc = tmp1 + tmp2 + tmp3;

o.Ksp = 10.^(log10Kspc)*exp(0.00019505*depth);%calcite mucci 1983

%Speciate
o.H = 10^-o.pH;
o.OH = o.Kw/o.H;
o.CO2 = o.DIC/(1+o.K1/o.H+o.K1*o.K2/(o.H^2));
o.HCO3 = o.DIC/(1+o.H/o.K1+o.K2/o.H);
o.CO3 = o.DIC/(1+o.H/o.K2+o.H^2/(o.K1*o.K2));
o.TA = o.CO2*(o.K1/o.H+2*o.K1*o.K2/o.H^2)+o.Kw/o.H-o.H;
o.chi = 1/(1+o.K2/o.H);
o.omega = o.Ca*o.CO3/o.Ksp;

%Forward k's
o.kf1 = o.prefactor*10^(329.85-110.541*log10(o.TK)-(17265.4/o.TK));     %Uchikawa and Zeebe (2012)

o.kf4 = o.prefactor*10^(13.635-2895/o.TK);                              %Uchikawa and Zeebe (2012)

%Backward k's
o.kb1 = o.kf1/o.K1;
o.kb4 = o.kf4*o.Kw/o.K1;
%Isotopologue concentrations - 14 total
o.EIC = o.CO3+o.HCO3;

%Seawater leak composition
o.DICsw = o.DIC;
o.Casw = o.Ca;
o.Srsw=o.Sr;
o.pHsw = o.pH;
o.Hsw = 10^-o.pHsw;
o.OHsw = o.Kw/o.Hsw;
o.CO2sw = o.DICsw/(1+o.K1/o.Hsw+o.K1*o.K2/(o.Hsw^2));
o.HCO3sw = o.DICsw/(1+o.Hsw/o.K1+o.K2/o.Hsw);
o.CO3sw = o.DICsw/(1+o.Hsw/o.K2+o.Hsw^2/(o.K1*o.K2));
o.Alksw = o.CO2sw*(o.K1/o.Hsw+2*o.K1*o.K2/(o.Hsw^2))+o.Kw/o.Hsw-o.Hsw;
o.EICsw = o.HCO3sw+o.CO3sw;

o.K_ABsw=o.K_AB;
o.K_MBsw=o.K_MB;
o.K_hco3sw=o.K_hco3;
%Cell composition
o.CO2cell = 13e-6;             %Chen et al. (2018);
end

function dy=carbonate(t,y,TK,kf1,kf4,kb1,kb4,K1,K2,Kw,Ksp,C266cell,C266sw,E2666sw,Casw,Srsw,Alksw,z,Dcell,tau,fCa,FAlk,k_A1,v_A1,v_M1,K_AB,K_MB,K_hco3)
dy=zeros(size(y));
DIC = y(1)+y(2);
KN = 0;
NT = 0;
TA = y(4);
a = 1;
b = KN+TA+K1;
c = TA*KN-KN*NT-Kw+TA*K1+KN*K1+K1*K2...
    -DIC*K1;
d = TA*KN*K1+TA*K1*K2-Kw*KN-KN*NT*K1-Kw*K1+KN*K1*K2...
    -DIC*KN*K1-2*DIC*K1*K2;
e = TA*KN*K1*K2-Kw*KN*K1-KN*NT*K1*K2-Kw*K1*K2...
    -DIC*KN*2*K1*K2;
f = -K1*K2*Kw*KN;
p = [a b c d e f];
r = roots(p);
H = max(real(r));
OH = Kw/H;
pH = -log10(H);
chi = 1/(1+K2/H);
omega = y(3)*(1-chi)*y(2)/Ksp;  %CO3 = y(4)(1-chi) 

A=y(3)*0.17;
M=y(5)*0.17;
phi1=-log10(K_hco3);
theta=10^(8.6-pH);
phi=10^(phi1-pH);

a0=1.7;
boltz=1.38065E-23;% Boltzmann constant

a_new=6.4e-10;
b_new=a_new/2;
h=3.1e-10;
d_new=27100;
alpha=1.41;
gamma=1.49e-10;
epsilon=10^(-19.9875);

B=omega*K_AB/A;


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


v_A2=v_A1;
v_A=v_A1+v_A2;

v_M2=v_M1;
v_M=v_M1+v_M2;


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
    
    Kp=rx/M*A;
    
    u_BA=k_B*B*P_A-v_B*P_B*(1-x);
    u_BM=k_BM*B*P_M-v_BM*P_B*x;
    
    u_net=u_A+u_M+u_BA+u_BM;
    
    Omega_A=2+(v_A*exp(2*epsilon/boltz/TK)/k_B/B+v_A*exp(2*epsilon/boltz/TK)*v_B/k_A/A/k_B/B)/(1-v_A*v_B/k_A/A/k_B/B);
    Omega_B=2+(v_B*exp(2*epsilon/boltz/TK)/k_A/A+v_A*exp(2*epsilon/boltz/TK)*v_B/k_A/A/k_B/B)/(1-v_A*v_B/k_A/A/k_B/B);
    
    C=2/u_net*(k_A*A/Omega_A+k_B*B/Omega_B);
    
    rho=((C^2+4*C)^(0.5)-C)/2;
    
    v_st=rho*u_net*b_new;
    sigma=abs(log(omega));
    beta_st=v_st/sigma;
    Delta_G=1/8*h*a_new*b_new*gamma*gamma/h/boltz/TK/sigma;
    I=beta_st*h/h/a_new/b_new*(4*h*sigma/3.14/h/a_new/b_new)^0.5*exp(-Delta_G/boltz/TK);
    Rnet = 1.137*h*(I*v_st^2)^(1/3)*d_new; % Calcite growth rate
  
end

JCaCO3=1e-3*Rnet;
x=rx/(1+rx);

dy(1) = -kf1*y(1)+kb1*y(2)*chi*H - kf4*y(1)*OH+kb4*y(2)*chi + Dcell/z*(C266cell-y(1)) + 1/tau*(C266sw-y(1));
dy(2) = kf1*y(1)-kb1*y(2)*chi*H + kf4*y(1)*OH-kb4*y(2)*chi + 1/tau*(E2666sw-y(2)) - 1/z*JCaCO3;
dy(3) = 1/tau*(Casw-y(3))+1/z*(0.5*fCa*FAlk/1000-JCaCO3*(1-x));
dy(4) = 1/tau*(Alksw-y(4))+1/z*(FAlk/1000-2*JCaCO3);
dy(5) = 1/tau*(Srsw-y(5)-+1/z*JCaCO3*x);
end


