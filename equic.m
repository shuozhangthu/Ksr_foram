% _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
%
% 	file: equic.m
%
%
%   Richard E. Zeebe* & Dieter A. Wolf-Gladrow**
%
%
%  *School of Ocean and Earth 
%   Science and Technology 
%   Department of Oceanography 
%   University of Hawaii at Manoa
%   1000 Pope Road, MSB 504
%   Honolulu, HI 96822
%   USA
%   email: zeebe@soest.hawaii.edu
%   http://www.soest.hawaii.edu/oceanography/faculty/zeebe.html
%
% 
% **Alfred Wegener Institute for       
%	Helmholtz Centre for Polar and Marine Research (AWI)	  
%	P.O. Box 12 01 61		       
%	D-27515 Bremerhaven		
%	Germany	      
%	e-mail: Dieter.Wolf-Gladrow@awi.de  	
%	http://www.awi-bremerhaven.de
%
% _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
%
% purpose: CO2 equilibrium constants. Compare book by
%
% Zeebe & Wolf-Gladrow, CO2 in Seawater: 
% Equilibrium, Kinetics, Isotopes.
%
%
% updates:
%
% 04/07/14 DWG & REZ: ocdflag for comparison with CO2SYS
%          KW P-coefficients water -> seawater
%          Kf only for ocdflag: ION (obsolete) -> iom0
% 03/26/14 REZ: equic.m: added p2f (pCO2 -> fCO2)
%          csys.m: output pCO2 and fCO2
%          R: 83.131 -> 83.14510. Kw: 148.96502 -> 148.9652 
%          (use old R to check old values)
% 03/02/10 REZ: boric acid pressure coefficient (10^3*a2)
%          corrected: 2.608 -> -2.608. Thanks to James Rae. 
%
% 15.02.01 density of seawater added
% 04.01.01 K1- and K2-Mehrbach as refit by Lueker added
% 04.01.01 phosphoric acid + pressure term added
% 03.01.01 k -> K
% 30.10.00 pressure eff. on Ksp* calcite and aragonite
% 06.12.99 Kf and pressure eff. on Kf
% 01.12.99 K_(calcite), K_(aragonite)
% 01.10.98 include pressure effect (Millero, 1995)
% 07.05.96 use constants summarized in Dickson and Goyet, 1994
%    
%
%
% MATLAB 5.3.1
% remarks:
%    units: mol/kg-soln
%    set phflag (0:Total), (1:Free)
% ===================================================================
tk = 273.15;           % [K] (for conversion [deg C] <-> [K])
T = TC + tk;           % TC [C]; T[K]
%Cl = S / 1.80655;      % Cl = chlorinity; S = salinity (per mille)
%cl3 = Cl.^(1/3);   
%ION = 0.00147 + 0.03592 * Cl + 0.000068 * Cl .* Cl;   % ionic strength
iom0 = 19.924*S/(1000.-1.005*S);
S_T = 0.14/96.062/1.80655*S;   % (mol/kg soln) total sulfate
                               %  Dickson and Goyet (1994) Ch.5 p.11

RGAS = 8.314510;        % J mol-1 deg-1 (perfect Gas)
R    = 83.14510;        % mol bar deg-1 
%R = 83.131;             % mol bar deg-1 
                        % conversion cm3 -> m3          *1.e-6
                        %            bar -> Pa = N m-2  *1.e+5
                        %                => *1.e-1 or *1/10                        

%------------------ Ks ----------------------------------------
%       Dickson and Goyet (1994), Chapter 5, p.13
%       (required for total2free)
%       Equilibrium constant for HSO4- = H+ + SO4--
%
%       K_S  = [H+]free [SO4--] / [HSO4-]
%       pH-scale: free scale !!!
%
%       the term log(1-0.001005*S) converts from
%       mol/kg H2O to mol/kg soln

tmp1 = -4276.1 ./ T + 141.328 -23.093*log(T);
tmp2 = +(-13856 ./ T + 324.57 - 47.986 * log(T)).*sqrt(iom0);
tmp3 = +(35474 ./ T - 771.54 + 114.723 * log(T)).*iom0;
tmp4 = -2698 ./ T .*sqrt(iom0).*iom0 + 1776 ./ T .*iom0 .*iom0;                                   

lnKs = tmp1 + tmp2 + tmp3 + tmp4 + log(1-0.001005*S);

Ks = exp(lnKs);

%------- total2free -----------------------------------------------
%
%       convert from pH_total ('total`) to pH ('free`):
%       pH_total = pH_free - log(1+S_T/K_S(s,tk))

total2free = 1.+S_T./Ks;

% --------------------- Kf  --------------------------------------------
%  Kf = [H+][F-]/[HF]  
%
%   (Dickson and Riley, 1979 in Dickson and Goyet, 
%   1994, Chapter 5, p. 14)
%   pH-scale: 'total'  
%
%   04/05/14 REZ: Kf only for ocdflag. ION (obsolete) -> iom0

tmp1 = 1590.2./T - 12.641 + 1.525.*sqrt(iom0);
tmp2 = log(1.-0.001005.*S) + log(1.+S_T./Ks);

if(ocdflag) % on free scale (omit last term)
 tmp2 = log(1.-0.001005.*S);
end;

lnKf = tmp1 + tmp2;
Kf = 0;

if phflag == 0;
        Kf  = exp(lnKf);
end;
if phflag == 1;
        lnKf = lnKf-log(total2free);
        Kf  = exp(lnKf);
end;

%------- sws2free -----------------------------------------------
%
%       convert from pH_sws ('seawater scale`) to pH ('free`):
%       pH_sws = pH_free - log(1+S_T/K_S(S,T)+F_T/K_F(S,T))

F_T = 7.e-5.*(S./35.);

if(ocdflag) % CO2SYS
 F_T = (0.000067./18.998).*(S./1.80655);
end;

sws2free   = (1.+S_T./Ks+F_T./Kf);

corr = sws2free./total2free;

% --------------------- Kwater -----------------------------------
%
%       Millero (1995)(in Dickson and Goyet (1994, Chapter 5, p.18))
%       $K_w$ in mol/kg-soln.
%       pH-scale: pH$_{total}$ ('total` scale).
%       
%       04/03/14 REZ: 148.96502 -> 148.9652                                                  

tmp1 = -13847.26./T + 148.9652 - 23.6521 .* log(T);
tmp2 = + (118.67./T - 5.977 + 1.0495.*log(T)).*sqrt(S) - 0.01615.*S;

if(ocdflag)
 % on sws (Dickson: +0.015) 
 tmp1 = tmp1 + 0.015;
end;

lnKw =  tmp1 + tmp2;

if phflag == 0;
        Kw  = exp(lnKw);
end;
if phflag == 1;
        lnKw = lnKw-log(total2free);
        Kw  = exp(lnKw);
end;

%---------------------- Kh (K Henry) ----------------------------
%
%               CO2(g) <-> CO2(aq.)
%               Kh      = [CO2]/ p CO2
%
%   Weiss (1974)   [mol/kg/atm]
%
%   Weiss + Murray & Riley (1971) data: 
%   T = 1-35 deg C, S = 0-38.
%                             
%

tmp = 9345.17 ./ T - 60.2409 + 23.3585 * log(T/100.);
nKhwe74 = tmp + S.*(0.023517-0.00023656*T+0.0047036e-4*T.*T);

%tmpKh1= 9050.69 ./ T - 58.0931 + 22.2940 * log(T/100.);
%nKhwe74l = tmpKh1 + S .* (0.027766-0.00025888*T+0.0050578e-4 * T .* T);
%Kh= exp(nKhwe74l);

Kh= exp(nKhwe74);

%---------------------- p2f (pCO2 -> fCO2) -----------------------
%               
%   convert from pCO2 to fCO2
%   Weiss (1974)
%

Pstd = 1.01325;
delC = (57.7 - 0.118.*T);
B = -1636.75 + 12.0408.*T - 0.0327957.*T.^2 + 3.16528.*0.00001.*T.^3;
p2f = exp((B + 2.*delC).*Pstd./(R*T));


% --------------------- K1 ---------------------------------------
%   first acidity constant:
%   [H^+] [HCO_3^-] / [CO2] = K_1
%
%   (Roy et al., 1993 in Dickson and Goyet, 1994, Chapter 5, p. 14)
%   pH-scale: 'total'. mol/kg-soln

tmp1 = 2.83655 - 2307.1266 ./ T - 1.5529413 .* log(T);
tmp2 =         - (0.20760841 + 4.0484 ./ T) .* sqrt(S);
tmp3 =         + 0.08468345 .* S - 0.00654208 .* S .* sqrt(S);   
tmp4 =         + log(1 - 0.001005 .* S);

lnK1roy = tmp1 + tmp2 + tmp3 + tmp4;
K1roy = 0;

if phflag == 0;
        K1roy  = exp(lnK1roy);
end;
if phflag == 1;
        lnK1roy = lnK1roy-log(total2free);
        K1roy   = exp(lnK1roy);
end;

% --------------------- K2 ----------------------------------------
%
%   second acidity constant:
%   [H^+] [CO_3^--] / [HCO_3^-] = K_2
%
%   (Roy et al., 1993 in Dickson and Goyet, 1994, Chapter 5, p. 15)
%   pH-scale: 'total'. mol/kg-soln

tmp1 = -9.226508 - 3351.6106 ./ T - 0.2005743 .* log(T);
tmp2 = (-0.106901773 - 23.9722 ./ T) .* sqrt(S);
tmp3 = 0.1130822 .* S - 0.00846934 .* S.^1.5 + log(1 - 0.001005 * S);

lnK2roy = tmp1 + tmp2 + tmp3;

if phflag == 0;
        K2roy  = exp(lnK2roy);
end;
if phflag == 1;
        lnK2roy = lnK2roy-log(total2free);
        K2roy   = exp(lnK2roy);
end;

% --------------------- K1 ---------------------------------------
%   first acidity constant:
%   [H^+] [HCO_3^-] / [H_2CO_3] = K_1
%
%   Mehrbach et al (1973) refit by Lueker et al. (2000).
%
%   Mehrbach data: T = 2-35 deg C, S = 19-43.
%
%   pH-scale: 'total'. mol/kg-soln

pK1mehr = 3633.86./T - 61.2172 + 9.6777.*log(T) - 0.011555.*S + 0.0001152.*S.*S;

if phflag == 0;
        K1mehr  = 10^(-pK1mehr);
end;
if phflag == 1;
	lnK1mehr = log(10^(-pK1mehr))-log(total2free);
        K1mehr   = exp(lnK1mehr);
end;

% --------------------- K2 ----------------------------------------
%
%   second acidity constant:
%   [H^+] [CO_3^--] / [HCO_3^-] = K_2
%
%   Mehrbach et al. (1973) refit by Lueker et al. (2000).
%
%   Mehrbach data: T = 2-35 deg C, S = 19-43.
%
%   pH-scale: 'total'. mol/kg-soln

pK2mehr = 471.78./T + 25.9290 - 3.16967.*log(T) - 0.01781.*S + 0.0001122.*S.*S;

if phflag == 0;
        K2mehr  = 10^(-pK2mehr);
end;
if phflag == 1;
	lnK2mehr = log(10^(-pK2mehr))-log(total2free);
        K2mehr   = exp(lnK2mehr);
end;

%----------- Roy or Mehrbach. default: Roy 

K1 = K1roy;
K2 = K2roy;

if (exist('k1k2flag'))

if k1k2flag == 0;
K1 = K1roy;
K2 = K2roy;
end;

if k1k2flag == 1;
K1 = K1mehr;
K2 = K2mehr;
end;

end;

% --------------------- Kb  --------------------------------------------
%  Kbor = [H+][B(OH)4-]/[B(OH)3]
%
%   (Dickson, 1990 in Dickson and Goyet, 1994, Chapter 5, p. 14)
%   pH-scale: 'total'. mol/kg-soln

tmp1 =  (-8966.90-2890.53*sqrt(S)-77.942*S+1.728*S.^(3./2.)-0.0996*S.*S);
tmp2 =   +148.0248+137.1942*sqrt(S)+1.62142*S;
tmp3 = +(-24.4344-25.085*sqrt(S)-0.2474*S).*log(T);

lnKb = tmp1 ./ T + tmp2 + tmp3 + 0.053105*sqrt(S).*T;

if phflag == 0;
        Kb  = exp(lnKb);
end;
if phflag == 1;
        lnKb = lnKb-log(total2free);
        Kb  = exp(lnKb);
end;

% --------------------- Phosphoric acid ---------------------
%
%
%   (DOE, 1994)  (Dickson and Goyet): pH_T, mol/(kg-soln)
%   Ch.5 p. 16
%
%

lnK1P = -4576.752 ./ T + 115.525 - 18.453*log(T) ...
        + (-106.736 ./ T + 0.69171) .* sqrt(S) ...
        + (-0.65643 ./ T - 0.01844) .* S;
lnK2P = -8814.715 ./ T + 172.0883 - 27.927 * log(T) ...
        + (-160.34 ./ T + 1.3566) .* sqrt(S) ...
        + (0.37335 ./ T - 0.05778) .* S;
lnK3P = -3070.75 ./ T - 18.141 ...
        + (17.27039 ./ T + 2.81197) .* sqrt(S) ...
        + (-44.99486 ./ T - 0.09984) .* S;

if(ocdflag)
 % on sws (Dickson: +0.015)
 lnK1P = lnK1P + 0.015;
 lnK2P = lnK2P + 0.015;
 lnK3P = lnK3P + 0.015;
end;
    
K1P = exp(lnK1P);
K2P = exp(lnK2P);
K3P = exp(lnK3P);

% --------------------- Silicic acid ---------------------------
%
%   (DOE, 1994)  (Dickson and Goyet): pH_T, mol/(kg-soln)
%   Ch.5 p. 17
%
%

%lnKSi = -8904.2 ./ T + 117.385 - 19.334*log(T) ...
%      + (3.5913-458.79 ./ T) .* sqrt(iom0) + (188.74 ./ T - 1.5998) .* iom0 ...
%      + (0.07871 - 12.1652 ./ T) .*iom0.^2 + log(1-0.001005*S);
%
%KSi = exp(lnKSi);

% --------------------- Kspc (calcite) ----------------------------
%
% apparent solubility product of calcite
%
%  Kspc = [Ca2+]T [CO32-]T
%
%  where $[]_T$ refers to the equilibrium total 
% (free + complexed) ion concentration.
%
%  Mucci 1983 mol/kg-soln

tmp1 = -171.9065-0.077993.*T+2839.319./T+71.595.*log10(T);
tmp2 = +(-0.77712+0.0028426.*T+178.34./T).*sqrt(S);
tmp3 = -0.07711.*S+0.0041249.*S.^1.5;
log10Kspc = tmp1 + tmp2 + tmp3;

Kspc = 10.^(log10Kspc);

% --------------------- Kspa (aragonite) ----------------------------
%
% apparent solubility product of aragonite
%
%  Kspa = [Ca2+]T [CO32-]T
%
%  where $[]_T$ refers to the equilibrium total 
% (free + complexed) ion concentration.
%
%  Mucci 1983 mol/kg-soln

tmp1 = -171.945-0.077993.*T+2903.293./T+71.595.*log10(T);
tmp2 = +(-0.068393+0.0017276.*T+88.135./T).*sqrt(S);
tmp3 = -0.10018.*S+0.0059415.*S.^1.5;
log10Kspa = tmp1 + tmp2 + tmp3;

Kspa = 10.^(log10Kspa);

%----------------------------------------------------
%
% Density of seawater as function of S,T,P.
%
% Millero et al. 1981, Gill, 1982.
%
%     
%                 
%
%----------------------------------------------------

%------------ Density of pure water

rhow = 999.842594 + 6.793952e-2*TC -9.095290e-3*TC^2 ...
            + 1.001685e-4*TC^3 -1.120083e-6*TC^4 + 6.536332e-9*TC^5;

%------------ Density of seawater at 1 atm, P=0

A =   8.24493e-1 - 4.0899e-3*TC + 7.6438e-5*TC^2 - 8.2467e-7*TC^3 ...
    + 5.3875e-9*TC^4;
    
B = -5.72466e-3 + 1.0227e-4*TC - 1.6546e-6*TC^2; 

C = 4.8314e-4;   

rho0 = rhow + A*S + B*S^(3/2) + C*S^2;

%-------------- Secant bulk modulus of pure water 
%
% The secant bulk modulus is the average change in pressure 
% divided by the total change in volume per unit of initial volume.

Ksbmw =   19652.21 + 148.4206*TC - 2.327105*TC^2 ...
	+ 1.360477e-2*TC^3 - 5.155288e-5*TC^4;

%-------------- Secant bulk modulus of seawater at 1 atm

Ksbm0 = Ksbmw ...
	+ S*( 54.6746 - 0.603459*TC + 1.09987e-2*TC^2 ...
			- 6.1670e-5*TC^3) ...
	+ S^(3/2)*( 7.944e-2 + 1.6483e-2*TC - 5.3009e-4*TC^2);


%-------------- Secant bulk modulus of seawater at S,T,P
	
Ksbm = Ksbm0 ...
	+ P*( 3.239908 + 1.43713e-3*TC + 1.16092e-4*TC^2 ...
		- 5.77905e-7*TC^3) ...
	+ P*S*( 2.2838e-3 - 1.0981e-5*TC - 1.6078e-6*TC^2) ...
 	+ P*S^(3/2)*1.91075e-4 ...
 	+ P*P*(8.50935e-5 - 6.12293e-6*TC + 5.2787e-8*TC^2) ...
 	+ P^2*S*(-9.9348e-7 + 2.0816e-8*TC + 9.1697e-10*TC^2);
 	
%------------- Density of seawater at S,T,P

rho = rho0/(1.-P/Ksbm);

if(ocdflag)
 % from total to sws at P=0
 K1 = K1*corr;
 K2 = K2*corr;
 Kb = Kb*corr;
 % Kw, KiP's already on sws for ocdflag = 1
end;

%---------------------- Pressure effect on K's (Millero, 95) ----------%
if(P > 0.0)

% index: K1 1, K2 2, Kb 3, Kw 4, Ks 5, Kf 6, Kspc 7, Kspa 8,
%        K1P 9, K2P 10, K3P 11

%----- note: there is an error in Table 9 of Millero, 1995.
%----- The coefficients -b0 and b1
%----- have to be multiplied by 1.e-3!

%----- there are some more errors! 
%----- the signs (+,-) of coefficients in Millero 95 do not
%----- agree with Millero 79

%----- Millero, 1995 Table 9: KW P-coefficients for water
%----- not seawater (cf. Millero, 1983). changed 04/05/14

a0 = -[25.5   15.82  29.48  20.02  18.03    9.78  48.76   46. ...
	14.51 23.12 26.57];
a1 =  [0.1271 -0.0219 0.1622 0.1119 0.0466 -0.0090 0.5304  0.5304 ...
	0.1211 0.1758 0.2020];
a2 =  [0.0     0.0   -2.608 -1.409 0.316  -0.942  0.0     0.0 ...
	-0.321 -2.647 -3.042]*1.e-3;
b0 = -[3.08   -1.13   2.84   5.13   4.53    3.91  11.76   11.76 ...
	2.67 5.15 4.08]*1.e-3;
b1 =  [0.0877 -0.1475 0.0    0.0794 0.09    0.054  0.3692  0.3692 ...
	0.0427 0.09 0.0714]*1.e-3;
b2 =  [0.0     0.0    0.0    0.0    0.0     0.0    0.0     0.0 ...
	0.0 0.0 0.0];

% Kw seawater/water
%a0(4) = -20.02; a1(4) = 0.1119; a2(4) = -1.4090e-3;
%a0(4) = -25.60; a1(4) = 0.2324; a2(4) = -3.6246e-3;

for ipc=1:length(a0);
  deltav(ipc)  =  a0(ipc) + a1(ipc).*TC + a2(ipc).*TC.*TC;
  deltak(ipc)  = (b0(ipc) + b1(ipc).*TC + b2(ipc).*TC.*TC);  
  lnkpok0(ipc) = -(deltav(ipc)./(R.*T)).*P + (0.5*deltak(ipc)./(R.*T)).*P.*P;
end;

K1 = K1*exp(lnkpok0(1));
K2 = K2*exp(lnkpok0(2));
Kb = Kb*exp(lnkpok0(3));
Kw = Kw*exp(lnkpok0(4));
Ks = Ks*exp(lnkpok0(5));
Kf = Kf*exp(lnkpok0(6));
Kspc = Kspc*exp(lnkpok0(7));
Kspa = Kspa*exp(lnkpok0(8));
K1P = K1P*exp(lnkpok0(9));
K2P = K2P*exp(lnkpok0(10));
K3P = K3P*exp(lnkpok0(11));

end; % Pressure effect

if(ocdflag)
 % from sws back to total at P=P
 corr = (1.+S_T./Ks+F_T./Kf)/(1.+S_T./Ks);
 K1 = K1/corr;
 K2 = K2/corr;
 Kb = Kb/corr;
 Kw = Kw/corr;
 K1P = K1P/corr;
 K2P = K2P/corr;
 K3P = K3P/corr;
end;

return;               