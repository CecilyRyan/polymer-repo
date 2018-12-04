%% Prediction of Binary Blend Morphology and CB localization
% this code predicts the morphology of a 2-phase polymer blend and the localization
% of a nano-particulate using the geometric mean equation, and contact
% angles of each polymer
% MIND YOUR UNITS
% Inputs: Contact angles of Polymer and liquid
% Outputs: interfacial tension, spreading coef
% EQNs are from CH3 of POLYMER INTERFACE AND ADHESION BY SOUHENG WO

clc 
clc 
clear
TI=190;%temperature the Blend is made at (C)
TM=20; %temperature Contact Angle is run at(C)
dt=TI-TM;

% Polymer(i)
namei='PHBV';
%Contace angels in H20 and DIIMETH
P1=(63.3+65.5+65.5+65.4+63.2)/5;
P2=(46.85+46.15+46.4+48.5+48.3)/5;
[Gi,Gid,Gip,Gim]= OwensWendt(P1, P2, TI, namei, TM); % SURFACE tension From owens wendt
dgdti=.06;%(mj/M^2C)



% Polymer(j)
namej='PLA';
Pj1=(63.15+63.25+64.7+64.4+64.5)/5;
Pj2=(65.7+59.55+65.8+59.15+59.9)/5;
[Gj,Gjd,Gjp,Gjm]= OwensWendt(Pj1, Pj2, TI, namej,TM); % SURFACE tension From owens wendt
dgdtj=.06;%(mj/M^2C)


%Filler CB 
namef='CB';
PF1=86.82;
PF2=59.13;
[Gcb,Gcbd,Gcbp,Gcbm]= OwensWendt(PF1, PF2, TI, namef,TM); % SURFACE tension From owens wendt
Gcb=98.1; % as measured at room temp 
Gcbd=84.1;% (200C, mJ/m^2) % From paper in polymers
Gcbp=3.2;% (200C, mJ/m^2)
dgdtcb=.06;%(mj/M^2C)
Gcbm=Gcb-(dgdtcb*(TI-TM));

% ---------------------Harmonic Mean EQN---------------------
% Inter_tension of i and j = Suf.ten(i) + Suf.ten(J)- Work of
% adhesion(dispersive -Work of adhesion(polar)

%Interfacial tension between polymers
Gij=Gi+Gj-(4*Gid*Gjd)/(Gid+Gjd)-(4*Gip*Gjp)/(Gip+Gjp);



%interfacial tension with CB
Gicb=Gi+Gcb-(4*Gid*Gcbd)/(Gid+Gcbd)-(4*Gip*Gcbp)/(Gip+Gcbp);
Gjcb=Gj+Gcb-(4*Gjd*Gcbd)/(Gjd+Gcbd)-(4*Gjp*Gcbp)/(Gjp+Gcbp);

%spreading coefficient 
Lij=Gj-Gi-Gij; %spreading Coef for PHBV on PLA 
Lji=Gi-Gj-Gij; % spreading coed for PLA on PHBV

% Contact angle Between two polymers

thetaij=acosd((Gj-Gij)/Gi);
thetaji=180-thetaij;

if Lij >0
    fprintf('%s will  spread at the interface of %s. \n',namei,namej)
elseif Lij <0 
    fprintf('%s will bead up at the interface of %s. \n',namei,namej)
end

if Lji >0
    fprintf('%s will spread at the interface of %s. \n\n',namej,namei)
elseif Lji <0 
    fprintf('%s will bead up at the interface of %s. \n\n',namej,namei)
end

Q=abs(Gi-Gj);
% Page 618 Introduction to physcial polymer science
if Gij > Q
    fprintf('Spreading coef will always be negative sujesting Immiscibility\n\n')
end

%Wetting coeff with CB
wij=(Gicb-Gjcb)/Gij; %PHBV-PLA
wji=(Gjcb-Gicb)/Gij;

if wij>1
    phasea= namej;
elseif wij <-1
    phasea=namei;
else 
    phasea='Interface';
end

if wji>1
    phase= namei;
elseif wji <-1
    phase=namej;
else 
    phase='Interface';
end

    
%Comparison Table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
Component={namei;namej;namef};% I J K CB
ST=[Gi;Gj;Gcb];
STmelt=[Gim;Gjm;Gcbm];
Dispersive=[Gid;Gjd;Gcbd];
Polar= [Gip;Gjp;Gcbp];
dgdt=[dgdti;dgdtj;dgdtcb];
TT= table( Component,ST, STmelt,Dispersive, Polar, dgdt);
disp(TT)
%
S='/';
IJ=[namei S namej];
IF=[namei S namef];
JF=[namej S namef];

JI=[namej S namei];
ComponentCouple={IJ;IF;JF};
InterfacialTension= [Gij;Gicb;Gjcb];
Blend={IJ;JI;0 };
SpreadingCoefficient=[Lij;Lji;0];
BlendCA=[thetaij; thetaji;0];
WettingCoefficient=[wij;wji;0];
Localization={phasea;phase;0};
T=table(ComponentCouple,InterfacialTension,Blend,SpreadingCoefficient,BlendCA, WettingCoefficient,Localization);
disp(T)