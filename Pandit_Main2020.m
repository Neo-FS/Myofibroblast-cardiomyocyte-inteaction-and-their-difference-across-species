% function [M1]=Pandit_Fusheng
% MODEL NAME: Pandit_Fusheng in rat
% SHORT DESCRIPTION: This model reproduces the action potential recorded experimentally
% for myocytes isolated from the adult rat left ventricle (Pandit et al. 2001) and
% active fibroblasts from Fusheng (Fusheng et al. 2021).
% Geometrical structure of left ventricular myocytes in rat are redesign
% baed on the the experimental measeuremnt (Wulfsohn et al, The Anatomical
% Record Part A: Discoveries in Molecular, Cellular, and Evolutionary
% Biology, 2004) and Luo report (Luo CH and  Rudy Y. A, Circulation Research. 1994)
% Becasue Pandit model cannot simulate the Ca2+ cycling of rat myocytes effectively,
% we replaced the Pandit model formulation for Ca2+ cycling with that of the
% ten Tusscher model (ten Tusscher et al. 2006).Additionally, Maximal Irel
% conductance, V_rel, is adjusted to 0.102 mM/ms, a transition rate of I_rel, k4
% is adjust to 0.005 ms-1,to agree well with the experimental data (Kaprielian et al, 1999)
% Maximum INaK current,I_NaK_max, and Scaling factor for INaCa , K_NaCa, are adjusted to
% satisfy the convergencein of Na+ and Ca2+ in epicardial myocytes and endocardial myocytes.
% Here 1 - Epicardical; 2 - Endocardial; 3 -TGF for Endocardial;

clc; clear all
t1=clock;
%------------------------------------------------------------------------
%                  Parameter for rat ventricular myocytes
%------------------------------------------------------------------------
Num_Fb = 0;  % the numeb of active fibroblast
Celltype = 2; % 1,Epicardial cell; 2,Endocardial cell
PD =[];PI =[];
%  Physical constant
R = 8314.5;
T = 295;
F = 96487;
Rtonf = (R*T)/F;
% Intracellular volumes
Cm = 0.132;% microF

% Geometrical structure of cardiomyocytes
Length = 0.11; width = 0.028; depth = width/3; % millimeter
Cell_volume = 1000*pi*Length*width*depth/4;%microliter
Vc = Cell_volume*0.68;%9.36e-3; %microliter-Myoplasm volume
Vsr = 0.0667*Vc; % microliter-Sarcoplasmic reticulum volume
Vss = 0.0033*Vc; % microliter-Submembrance volume
inverseVcF = 1/(Vc*F);
inverseVcF2 = 1/(2*Vc*F);
inversevssF2 = 1/(2*Vss*F);

% Extracellular concentrations
Ca_o = 1.2;
Na_o=140;
K_o = 5.4;

if  Celltype == 1
    G_Na=0.8;   % pA/pF
    G_to=0.035;  % pA/pF
    a_to = 0.886;
    b_to = 0.114;
    K_NaCa=0.000009984*0.7525;
    I_NaK_max=0.08*1.12;
    alpha_TGF = 0;
elseif  Celltype == 2
    G_Na=0.8*1.3;   % pA/pF
    G_to=0.035*0.5;   % pA/pF
    a_to = 0.583;
    b_to = 0.417;
    K_NaCa=0.000009984*0.94;
    I_NaK_max=0.08*1.24;
    alpha_TGF = 0;
elseif Celltype ==3
    G_Na=0.8*1.3*1.39045;%TGF level change the INa current;
    G_to=0.035*0.57404;   %TGF level change the Ito current;
    a_to = 0.886;
    b_to = 0.114;
    K_NaCa=0.000009984*1.00;
    I_NaK_max=0.08*1.24;
    alpha_TGF = 0.3;
else
    error('Error:Incorrect input of parameters')
end
% Parameters for I_NaCa
K_m_K=1.5;
K_m_Na=10;
I_pCa_max=0.004;
d_NaCa=0.0001;
gamma_NaCa=0.5;

G_CaL=0.031;
E_CaL=65;
tau_Ca_inact=9;

%   Parameters for I_ss
G_ss=0.007;
%  Parameters for I_K1
G_K1=0.024;
%  Parameters for I_f
G_f=0.00145;
f_Na=0.2;
%   Parameters for Background_currents
G_bNa=0.00008015;
G_bCa=0.0000324;
G_bK=0.000138;
% Parameters for I_pCa
G_pCa=0.01238;
K_pCa=0.0005;

% Calcium buffering dynamics
Bufc=0.20;
Kbufc=0.001;
Bufsr=10.;
Kbufsr=0.3;
Bufss=0.40;
Kbufss=2.5e-4;
Vmaxup=6.375e-3;
Kup = 2.5e-4;
Vrel = 0.102;%%40.8;%%
k1_= 0.15;
k2_= 0.045;
k3 = 0.060;
k4 = 0.005;%%1.5e-5;%%
EC=1.5;
maxsr = 2.5;
minsr = 1.;
Vleak = 3.6e-4;
Vxfer = 3.8e-3;

k_plus_htrpn=0.00237e3;
k_minus_htrpn=3.2e-5;
k_plus_ltrpn=0.0327e3;
k_minus_ltrpn=0.0196;
K_Ca =0.5994;

LTRPN_tot = 0.07;
HTRPN_tot = 0.14;
LTRPN_Ca=0.011163228644803;
HTRPN_Ca=0.130189159299828;

% Initial intracellular  concentrations of cardiomyocyte
Na_i=10.7347;
K_i =139.2706;
Ca_i = [0.000116520413762464];
Ca_ss = [0.000166615985336238];
Ca_SR = [3.23401918550965];
% Initial value of ion channge gate in myocytes
V = [-80.3747936704937];
m_Na = [0.00431095351450453];
h_Na = [0.669108788377134];
j_Na = [0.668969383741482];
r_ss = [0.00293833925906765];
s_ss = [0.310015061307104];
y_f = [0.00347188180148011];
r_to = [0.00221590376703238];
s_to = [0.704802198192475];
s_slow = [0.485453218122925];
d_Ca = [2.22676681555124e-06];
f_1Ca = [0.999951782141957];
f_2Ca = [0.999951782159783];
Ca_inact = [0.999888336732021];
sRR = [0.992667363885394];

%--------------------------------------------------------------------------
%               PARAMETERS FOR INITIAL  Contraction
%--------------------------------------------------------------------------
TRPNCa = [0.00115400000000000];
N0=[0.999238085847744];
N1=[1.75743609919447e-05];
P0=[7.22780384324895e-05];
P1=[0.000195319631858015];
P2=[0.000209241625239671];
P3=0.523108e-3;
SL=1.9;
SL_0=1.9;
k_PN=0.045;
f_XB=0.10;
g_minxb=0.14;

%--------------------------------------------------------------------------
%             Constants of myofibroblast Parameter
%--------------------------------------------------------------------------
% Extracelluer concentrations of myofibroblast
fNa_o=140.0;%
fK_o=5.4;%
% Intacelluer  concentrations of myofibroblast
fNa_i=10; % 8.5547 ;mM;
fK_i= 129.4349; %mM;

% fI_Na
GNa_f = 4.2;

%  fI_b
fG_bNa=4.0e-3;
fG_bK=4.96e-3;

% fI_MGC
G_fMGC = 0.018;%pS/pF
detaCL = 0;

%  fI_NaK
fV_rev=-150.0;
fB=-200;
fK_mk=1.0;%mM3.0
fK_mNa=11.0;%mM
fI_NaKoo=1.75;

%   I_K1
fG_K1=0.178;
a_K1=0.94;
b_K1=1.26;
fC_m = 18;%  pF;
G_gap = 7;%  nS

%  Initial Values of fibroblast model for [K_o]=5.4mM
P_to = 3.540000000000000e-04;
kvo=30e-3;%
k_vo=2e-3;%
zv=1.28;%
z_v=-1.53;%
ko=77e-3;%
k_o=18e-3;%

m_f = 0.166551493609004;
h_f = 0.021549092624088;
j_f = 0.021549092624086;
C_0to = 0.002744277836910;
C_1to = 4.100094085704080e-04;
C_2to = 2.297157828581396e-05;
C_3to = 5.720122491618227e-07;
C_4to = 5.341350664715740e-09;
O_to = 2.284911117674759e-08;
V_f = -54.238032147101580;

%--------------------------------------------------------------------------
%             PARAMETER FOR SIMULATION DURATION
%--------------------------------------------------------------------------
% for stims2 = 400:50:500;
dt=0.005;
xns = 5; % The number of cycle
disp(['The number of cycle: ',num2str(xns)]);
sdur = 1.0;
stimstrength= -6;  %  stimstrength
stims1 = 1000;
stims2 = 1000;
tbegin=10;
xnstims1 = stims1+tbegin+sdur;%+tbegin;     ???
xnstims2 =stims1*(xns-1)+tbegin+stims2+sdur;%tbegin;    ????
endtime=xnstims2+stims2;  %duration of the simulation
nswitch = 0;
ncounts1=1;
napdd=0;
apdtime = zeros(10,1);

%    time of the simulation
time=0; step=0;
sst=[];ssv=[];ssa=[];ssb=[];ssc=[];ssd=[];sse=[];ssf=[];ssg=[];
APV = [];
while time<=endtime
    %--------------------------------------------------------------------------
    %            Simulation protocols
    %--------------------------------------------------------------------------
    if (time >=tbegin&&time<=tbegin+sdur)
        I_stim = stimstrength;
    elseif (time >=xnstims1-sdur&&time<=xnstims1)
        I_stim = stimstrength;
        nswitch=1;
    else
        I_stim = 0;
    end
    
    if (time>xnstims1&& ncounts1 < xns-1 && nswitch==1)
        ncounts1=ncounts1+1;
        xnstims1=ncounts1*stims1+tbegin+sdur;
        nswitch=0;
    end
    
    if (time >=xnstims2-sdur&&time<=xnstims2)
        I_stim = stimstrength;
        skip=1;
        ncount=0;
    end
    
    %   The inversion of potential
    E_Na = Rtonf*log(Na_o/Na_i);
    E_Ca = 0.5*Rtonf*log(Ca_o/Ca_i);
    E_K = Rtonf*log(K_o/K_i);
    
    %  Fast Na+ current I_Na
    m_inf=(1/(1+exp((V+(45 ))/((-1)*(6.5 )))));
    tau_m=((1.36 )/((((0.32 )*(V+(47.13 )))/(1-exp(((-1)*0.1)*(V+(47.13 )))))+(0.08*exp(((-1)*V)/(11 )))));
    h_inf=(1/(1+exp((V+(76.1 ))/(6.07 ))));
    j_inf=h_inf;
    if V>=-40
        tau_h = 0.4537*(1+exp(-(V+10.66)/11.1));
        tau_j = 11.63 *(1+exp(-0.1*(V+32)))/exp(-2.535E-7*V);
    else
        tau_h = 3.49/(0.135*exp(-(V+80)/6.8)+3.56*exp(0.079*V)+310000*exp(0.35*V));
        tau_j =3.49/((((V+37.78)/(1+exp(0.311*(V+79.23))))*...
            (-127140*exp(0.2444*V)-3.474E-5*exp(-0.04391*V)))...
            +((0.1212*exp(-0.01052*V))/(1+exp(-0.1378*(V+40.14)))));
    end
    m_Na = m_inf- (m_inf-m_Na)*exp(-dt/tau_m);
    h_Na = h_inf-(h_inf-h_Na)*exp(-dt/tau_h);
    j_Na = j_inf-(j_inf-j_Na)*exp(-dt/tau_j);
    I_Na = G_Na*(m_Na^3)*h_Na*j_Na*(V-E_Na);
    
    %   Background_currents
    I_bNa = G_bNa*(V-E_Na);
    I_bCa = G_bCa*(V-E_Ca);
    I_bK = G_bK*(V-E_K);
    I_b = I_bNa+I_bCa+I_bK;
    
    %   Transient outward K currents I_to
    r_inf = 1/(1+exp((V+10.6)/-11.42 ));
    tau_r = 1000/((45.16*exp(0.03577*(V+50)))+(98.9*exp(-0.1*(V+38 ))));
    s_inf = 1/(1+exp((V+45.3)/(6.8841)));
    tau_s = (350*exp(-(((V+70)/15)^2)))+35;
    s_slow_inf=s_inf;
    if   Celltype==1
        tau_s = 350*exp(-(V+70)/15 )^2+35;
        tau_s_slow = 3700*exp(-((V+70)/30)^2)+35;
    elseif Celltype==2|Celltype == 3
        tau_s= 550*exp(-(((V+70 )/25 )^2))+49;
        tau_s_slow = 3300*exp(-((V+70)/30)^2)+49;
    else
        error('Error:Incorrect input of parameters')
    end
    r_to=r_inf-(r_inf-r_to)*exp(-dt/tau_r);
    s_to=s_inf-(s_inf-s_to)*exp(-dt/tau_s);
    s_slow = s_slow_inf-(s_slow_inf-s_slow)*exp(-dt/tau_s_slow);
    I_to = G_to*r_to*(a_to*s_to+b_to*s_slow)*(V-E_K);
    
    %   Steady State Outward K currents I_sus
    r_ss_inf = 1/(1+exp((V+(11.5))/(-11.82)));
    tau_r_ss = 10000/((45.16*exp(0.03577*(V+50)))+(98.9*exp(-0.1*(V+38))));
    s_ss_inf = 1/(1+exp((V+87.5)/(10.3)));
    tau_s_ss = 2100;
    r_ss = r_ss_inf-(r_ss_inf-r_ss)*exp(-dt/tau_r_ss);
    s_ss = s_ss_inf-(s_ss_inf-s_ss)*exp(-dt/tau_s_ss);
    I_ss = G_ss*r_ss*s_ss*(V-E_K);
    
    %  Hyperpolarisation Active I_f current
    y_inf = (1/(1+exp((V+138.6)/(10.48))));
    tau_y = 1000/((0.11885*exp((V+80)/(28.37)))+(0.5623*exp((V+80)/(-14.19))));
    y_f = y_inf-(y_inf-y_f)*exp(-dt/tau_y);
    
    f_K = (1-f_Na);
    I_fNa = (((G_f*y_f )*f_Na)*(V-E_Na));
    I_fK = (((G_f*y_f )*f_K)*(V-E_K));
    I_f = I_fNa+ I_fK;
    
    %   Inward_rectifer I_K1 current
    I_K1 = (48/(exp((V+37)/25)+exp(-(V+37)/25))+10)*(0.001/(1+exp((V-E_K-76.77)/-17)))+...
        G_K1*(V-E_K-1.73)/((1+exp(1.613*(V-E_K-1.73)/Rtonf))*(1+exp(-(K_o-0.9988)/0.124)));
    
    % sodium_potassium_pump
    sigma=((exp(Na_o/(67.3))-1)/7);
    I_NaK=((((((I_NaK_max*1)/((1+(0.1245*exp(((((-1)*0.1)*V)*F)/(R*T))))+((0.0365*sigma)*exp((((-1)*V)*F)/(R*T)))))*K_o)/(K_o+K_m_K))*1)/(1+((K_m_Na/Na_i)^1.5)));
    
    % Na_Ca_ion_exchanger_current
    I_NaCa=((K_NaCa*((((Na_i^3)*Ca_o)*exp((0.03743*V)*gamma_NaCa))-(((Na_o^3)*Ca_i)*exp((0.03743*V)*(gamma_NaCa-1)))))/(1+(d_NaCa*((Ca_i*(Na_o^3))+(Ca_o*(Na_i^3))))));
    
    %    Sarcolemmal_calcium_pump_current
    I_pCa=G_pCa*Ca_i/(K_pCa+Ca_i);
    
    %   L_type_Ca_channel  current
    d_inf = 1/(1+exp((V+15.3)/-5));
    tau_d = 3.05*exp(-0.0045*((V+7)^2))+1.05*exp(-0.002*((V-18)^2))+0.25;
    f11_inf = 1/(1+exp((V+(26.7))/5.4));
    tau_f_1Ca = 105*exp(-(((V+45)/12)^2))+40/(1+exp((-V+25)/25))+15/(1+exp((V+75)/25))+1.7;
    tau_f_2Ca = 41*exp(-(((V+47)/12)^2))+80/(1+exp((V+55)/(-5)))+15/(1+exp((V+75)/25))+1.7;
    Ca_inact_inf= 1/(1+(Ca_ss/1.5));
    d_Ca = d_inf-(d_inf-d_Ca)*exp(-dt/tau_d);
    f_1Ca= f11_inf- (f11_inf-f_1Ca)*exp(-dt/tau_f_1Ca);
    f_2Ca= f11_inf-(f11_inf-f_2Ca)*exp(-dt/tau_f_2Ca);
    Ca_inact = Ca_inact_inf-(Ca_inact_inf-Ca_inact)*exp(-dt/tau_Ca_inact);
    I_CaL= G_CaL*d_Ca*((0.9+Ca_inact/10)*f_1Ca+(0.1-Ca_inact/10)*f_2Ca)*(V-E_CaL);
    
    %    update concentrationsk
    Na_i = Na_i-dt*(I_Na+I_bNa+I_NaCa*3+I_NaK*3+I_fNa)*inverseVcF;
    K_i = K_i-dt*(I_stim+I_ss+I_bK+I_to+I_K1+I_fK-2*I_NaK)*inverseVcF;
    
    dLTRPN_Ca = k_plus_ltrpn*Ca_i*(LTRPN_tot-LTRPN_Ca)-k_minus_ltrpn*LTRPN_Ca;
    LTRPN_Ca = LTRPN_Ca+dt*dLTRPN_Ca;
    dHTRPN_Ca = k_plus_htrpn*Ca_i*(HTRPN_tot-HTRPN_Ca)-k_minus_htrpn*HTRPN_Ca;
    HTRPN_Ca = HTRPN_Ca+dt*dHTRPN_Ca;
    J_trpn = dLTRPN_Ca+dHTRPN_Ca;
    
    kCa_SR=maxsr-((maxsr-minsr)/(1+(EC/Ca_SR)*(EC/Ca_SR)));
    k1=k1_/kCa_SR;k2=k2_*kCa_SR;
    dRR=k4*(1-sRR)-k2*Ca_ss*sRR;
    sRR=sRR+dt*dRR;
    sOO=k1*Ca_ss*Ca_ss*sRR/(k3+k1*Ca_ss*Ca_ss);
    Irel=Vrel*sOO*(Ca_SR-Ca_ss)*(1-alpha_TGF);
    Ileak=Vleak*(Ca_SR-Ca_i);
    Iup=Vmaxup/(1.+((Kup*Kup)/(Ca_i*Ca_i)));
    Ixfer=Vxfer*(Ca_ss-Ca_i);
    CaCSQN=Bufsr*Ca_SR/(Ca_SR+Kbufsr);
    dCa_SR=dt*(Iup-Irel-Ileak);
    bjsr=Bufsr-CaCSQN-dCa_SR-Ca_SR+Kbufsr;
    cjsr=Kbufsr*(CaCSQN+dCa_SR+Ca_SR);
    Ca_SR=(sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;
    Ca_ssBuf=Bufss*Ca_ss/(Ca_ss+Kbufss);
    dCa_ss=dt*(-Ixfer*(Vc/Vss)+Irel*(Vsr/Vss)+(-I_CaL*inversevssF2));
    bcss=Bufss-Ca_ssBuf-dCa_ss-Ca_ss+Kbufss;
    ccss=Kbufss*(Ca_ssBuf+dCa_ss+Ca_ss);
    Ca_ss=(sqrt(bcss*bcss+4*ccss)-bcss)/2;
    CaBuf=Bufc*Ca_i/(Ca_i+Kbufc);
    dCa_i=dt*((-(I_bCa+I_pCa-2*I_NaCa)*inverseVcF2)-(Iup+J_trpn-Ileak)*(Vsr/Vc)+Ixfer);
    bc=Bufc-CaBuf-dCa_i-Ca_i+Kbufc;
    cc=Kbufc*(CaBuf+dCa_i+Ca_i);
    Ca_i=(sqrt(bc*bc+4*cc)-bc)/2;
    
    %%   Contraction
    f_01=3*f_XB;
    f_12=10*f_XB;
    f_23=7*f_XB;
    SL_norm=(SL-1.3)/(2.3-1.3);
    g_xbSL=g_minxb*(2-(SL_norm)^1.6);
    g_10SL=g_xbSL;
    g_21SL=2*g_xbSL;
    g_32SL=3*g_xbSL;
    
    Ntm=5+3*SL_norm;
    K_half=1/(1+K_Ca/(1.5-SL_norm));
    k_NP=k_PN*(LTRPN_Ca/(LTRPN_tot*K_half))^Ntm;
    
    N0=N0+dt*(k_PN*P0+g_10SL*N1-k_NP*N0);
    N1=N1+dt*(k_PN*P1-(k_NP+g_10SL)*N1);
    P0=P0+dt*(g_10SL*P1+k_NP*N0-(k_PN+f_01)*P0);
    P1=P1+dt*(g_21SL*P2+f_01*P0+k_NP*N1-(f_12+g_10SL+k_PN)*P1);
    P2=P2+dt*(g_32SL*P3+f_12*P1-(f_23+g_21SL)*P2);
    P3=P3+dt*(f_23*P2-g_32SL*P3);
    
    Sigma=g_minxb*2*g_minxb*3*g_minxb+f_01*2*g_minxb*3*g_minxb+f_01*f_12*3*g_minxb+f_01*f_12*f_23;
    P1_max=(f_01*(2*g_minxb)*(3*g_minxb))/(Sigma);
    P2_max=(f_01*f_12*3*g_minxb)/(Sigma);
    P3_max=(f_01*f_12*f_23)/(Sigma);
    F_max=P1_max+2*P2_max+3*P3_max;
    
    F_contrn=-(P1+N1+2*P2+3*P3)/(F_max);
    SL=0.8*F_contrn+SL_0;
    F_contr=-120*F_contrn;
    
    %--------------------------------------------------------------------------
    %           Compute Myofibbroblasts
    %--------------------------------------------------------------------------
    fE_Na = Rtonf*log(fNa_o/fNa_i);
    fE_K=Rtonf*log(fK_o/fK_i);
    
    fm_inf = 1.0/(1.0+exp(-(V_f+42)/7.6));
    fh_inf=1.0/(1.0+exp((V_f+84)/7.8));
    tm=1.36/((0.32*(V_f+47))/(1.0-exp(-0.1*(V_f+47)))+0.08*exp(-V_f/11));
    th = 1.0/(1.03e-5*exp(-(V_f+4.5)/7.5)+0.623*exp((V_f+0.5)/73.5));
    tj =1.0/(0.00091*exp(-(V_f+61)/9.6)+0.057*exp((V_f+108)/334));
    m_f = fm_inf -(fm_inf -m_f)*exp(-dt/tm);
    h_f = fh_inf-(fh_inf-h_f)*exp(-dt/th);
    j_f = fh_inf-(fh_inf-j_f)*exp(-dt/tj);
    fI_Na = GNa_f*m_f*h_f*j_f*(V_f-fE_Na);
    
    %  I_K1
    O_K1=1/(a_K1+exp(b_K1*(V_f-fE_K)/Rtonf));
    fI_K1=fG_K1*O_K1*sqrt(fK_o*0.001)*(V_f-fE_K);
    
    kv=kvo*exp(V_f*zv/Rtonf);
    k_v=k_vo*exp(V_f*z_v/Rtonf);
    if C_0to<0
        C_0to=0;
    elseif C_0to>1
        C_0to=1;
    else
        dC_0to=(dt)*(k_v*C_1to-4*kv*C_0to);
        C_0to=C_0to+dC_0to;
    end
    if C_1to<0
        C_1to=0;
    elseif C_1to>1
        C_1to=1;
    else
        dC_1to=(dt)*(2*k_v*C_2to+4*kv*C_0to-(3*kv+k_v)*C_1to);
        C_1to=C_1to+dC_1to;
    end
    if C_2to<0
        C_2to=0;
    elseif C_2to>1
        C_2to=1;
    else
        dC_2to=(dt)*(3*k_v*C_3to+3*kv*C_1to-(2*kv+2*k_v)*C_2to);
        C_2to=C_2to+dC_2to;
    end
    if C_3to<=0
        C_3to=0;
    elseif C_3to>=1
        C_3to=1;
    else
        dC_3to=(dt)*(4*k_v*C_4to+2*kv*C_2to-(kv+3*k_v)*C_3to);
        C_3to=C_3to+dC_3to;
    end
    if C_4to<=0
        C_4to=0;
    elseif C_4to>=1
        C_4to=1;
    else
        dC_4to=(dt)*(k_o*O_to+kv*C_3to-(ko+4*k_v)*C_4to);
        C_4to=C_4to+dC_4to;
    end
    if O_to<0
        O_to=0;
    elseif O_to>1
        O_to=1;
    else
        dO_to=(dt)*(ko*C_4to-k_o*O_to);
        O_to=O_to+dO_to;
    end
    fI_to=P_to*O_to*(V_f*F/Rtonf)*(fK_i-fK_o*exp(-V_f/Rtonf))/...
        (1-exp(-V_f/Rtonf));
    %  fI_b & fI_NaK
    fI_b=fG_bK*(V_f-fE_K)+fG_bNa*(V_f-fE_Na);
    fI_NaK=fI_NaKoo*(fK_o/(fK_o+fK_mk))*((V_f-fV_rev)/(V_f-fB))*...
        ((fNa_i/(fNa_i+fK_mNa))^1.5);
    % fI_MGC
    if detaCL == 0
        O_fMGC = 0;
    elseif detaCL > 4.02
        O_fMGC = 0;
    else
        O_fMGC = 0.0415*(detaCL^2.286);
    end
    fI_MGC = G_fMGC*O_fMGC*V_f;
    %--------------------------------------------------------------------------
    %                              Write the PARAMEITER
    %--------------------------------------------------------------------------
    I_inter_fibro=(V_f-V)*G_gap;
    I_inter_myo= -Num_Fb *I_inter_fibro;
    fI_tot = fI_Na+fI_K1+fI_to+fI_b+fI_NaK+ fI_MGC;
    I_tot = (I_Na+I_CaL+I_to+I_ss+I_f+I_K1+I_b+I_NaK+I_NaCa+I_pCa+I_stim)/Cm;
    if Num_Fb == 0
        V_f=V_f-dt*fI_tot;
    else
        V_f=V_f-dt*(fI_tot+I_inter_fibro/fC_m);
    end
    svol=V-dt* (I_tot+I_inter_myo/132);
    
    
    if(time > (stims1*(xns-2)))&& (time < (stims1*(xns-1)+tbegin))
        APV = [APV; V];
    end
    
    if(time > (stims1*(xns-1)+tbegin))&& (time < (stims1*(xns-1)+sdur+tbegin))
        %  AP50=max(APV)-(max(APV)-min(APV))*0.5;
        AP90=max(APV)-(max(APV)-min(APV))*0.9;
    end
    
    if(time > (stims1*(xns-1)+sdur+tbegin))
        if ((svol-AP90)*(V-AP90)<0)
            napdd=napdd+1;
            apdtime(napdd) = time;
        end
        if(napdd==3)
            DI=apdtime(2)-apdtime(1);
            APD=apdtime(3)-apdtime(2);
        end
    end
    V = svol;
    
    if mod(step,50)==0 % To improve the speed of computer
        sst = [sst;time];
        ssv = [ssv; V];
        ssa = [ssa; V_f];
        ssb = [ssb; F_contr];
        ssc = [ssc; SL];
        ssd = [ssd; Ca_i*1000]; % To change the units
        sse = [sse; Na_i];
        ssf = [ssf; K_i];
    end
    step=step+1;
    time=time+dt;
end
AA = [apdtime(1) AP90; apdtime(2) AP90];
BB = [apdtime(2) AP90; apdtime(3) AP90];
if(napdd==3)
    PD =[PD;APD];
    PI =[PI; DI];
end
% end
figure
subplot(2,2,1)
plot(sst,ssv,'linewidth',2)
hold on
plot(BB(:,1),BB(:,2),'r','linewidth',2)
xlabel('time (ms)','FontSize',15)
ylabel('Cardiomtyocyte voltage (mV)','FontSize',15)
hold on
plot(AA(:,1),AA(:,2),'g','linewidth',2)
hold on
legend('Action potential',['APD = ',num2str(PD)],['DI = ',num2str(DI)])

subplot(2,2,2)
plot(sst,ssa,'linewidth',2)
xlabel('time (ms)','FontSize',15)
ylabel('Myofibroblast voltage (mV)','FontSize',15)

subplot(2,2,3)
plot(sst,ssb,'linewidth',2)
xlabel('time (ms)','FontSize',15)
ylabel('Active tension (kPa)','FontSize',15)

subplot(2,2,4)
plot(sst,ssc,'linewidth',2)
xlabel('time (ms)','FontSize',15)
ylabel('Active deformation (kPa)','FontSize',15)

figure
subplot(2,2,1)
plot(sst,ssd,'linewidth',2)
xlabel('time (ms)','FontSize',15)
ylabel('Cardiomtyocyte Ca_i (nM)','FontSize',15)

subplot(2,2,2)
plot(sst,sse,'linewidth',2)
xlabel('time (ms)','FontSize',15)
ylabel('Cardiomtyocyte Na_i (uM)','FontSize',15)

subplot(2,2,3)
plot(sst,ssf,'linewidth',2)
xlabel('time (ms)','FontSize',15)
ylabel('Cardiomtyocyte K_i (uM)','FontSize',15)

M2=[sst,ssv,ssa,ssb,ssc,ssd,sse,ssf];
dlmwrite('Pandit_Myofibroblast.txt',M2,'delimiter','\t','precision',6)
t2 = clock;
disp(['total time: ',num2str(etime(t2,t1))]);

