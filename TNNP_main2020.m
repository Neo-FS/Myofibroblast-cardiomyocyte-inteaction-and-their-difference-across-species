% MODEL NAME:TNNP_MacCannell in humans
% SHORT DESCRIPTION: This model reproduces the action potential recorded experimentally
% for myocytes isolated from the adult human left ventricle (Ten Tusscher et al. 2006) and
% active fibroblasts from Fusheng (Fusheng et al. 2021).
% Maximum INaK current,I_NaK_max, and Scaling factor for INaCa , K_NaCa, are adjusted to
% satisfy the convergencein of Na+ and Ca2+ in epicardial myocytes, M myocytes and endocardial myocytes.

% Here 1 - Epicardical; 2 - Endocardial; 3 - M cardial£» 4 -TGF for Endocardial;

clc;clear all
t1=clock;
%------------------------------------------------------------------------
%                  Parameter for human ventricular myocytes
%------------------------------------------------------------------------
Num_Fb = 0;  % the numeb of active fibroblast
Celltype = 4; % 1,Epicardial cell; 2,Endocardial cell;  3, Middle
%  Physical constant
R=8314.5;  % millijoule_per_mole_kelvin (in membrane) - Ideal gas constant
T=310;     % kelvin (in membrane) - Absolute tempehumanure
F=96487;   % coulomb_per_mole (in membrane) - Faraday constant
Rtonf = (R*T)/F;

%  Cellular capacitance
C_m=0.185; % microF (in membrane) - human myocyte membrane capacitance

% A=load('TNNP2004data.txt');
%------------------------------------------------------------------------
%                       the ventricular myocytes
%------------------------------------------------------------------------
%External concentrations
v=-86.2;
K_o=5.4;
Ca_o=2.0;
Na_o=140.0;
Ca_i = 9.558209159100939e-05;
CaSS = 2.063017063722722e-04;
CaSR = 3.181497637568659;%2.828640148792239;
Na_i= 7.777862482160443;%   mmol/L;
K_i=1.376172340379077e+02;%   mmol/L;

%% Intracellular
Cm=0.185; %nF
R=8314.472;
F=96485.3415;
T=295.0;
Rtonf=(R*T)/F;

Vc=0.016404;
Vsr=0.001094;
Vss=0.00005468;
inverseVcF2=Cm/(2*Vc*F);
inverseVcF=Cm/(Vc*F);
inversevssF2=Cm/(2*Vss*F);

k_plus_htrpn=0.00237e3;
k_minus_htrpn=3.2e-5;
k_plus_ltrpn=0.0327e3;
k_minus_ltrpn=0.0196;
% K_Ca=k_minus_ltrpn/k_plus_ltrpn;
K_Ca =0.5994;
LTRPN_tot = 0.07;
HTRPN_tot = 0.14;
LTRPN_Ca=0.0112684;
HTRPN_Ca=0.12529;

% Calcium buffering dynamics
Bufc=0.17;
Kbufc=0.001;
Bufsr=10.;
Kbufsr=0.3;
Bufss=0.4;
Kbufss=0.00025;

% Intracellular calcium flux dynamics
Vmaxup=0.006375;
Kup=0.00025;
Vrel=0.102;%%40.8;%%
k1_=0.15;
k2_=0.045;
k3=0.060;
k4=0.005;%%0.000015;%%
EC=1.5;
maxsr=2.5;
minsr=1.;
Vleak=0.00036;
Vxfer=0.0038;

% Parameters for I_Kr  currents
G_Kr = 0.153;
%%  Parameters for and I_to Iks
alpha=2.5;
p_KNa=0.03;
if Celltype == 1
    G_Ks = 0.392;  % nanoS per picoF;
    G_to = 0.294;  % nanoS per picoF;
    knaca=1030.0;
    P_NaK=2.724*1.20;
    alpha_TGF = 0;
elseif Celltype == 2
    G_Ks = 0.392;
    G_to = 0.073;
    knaca=980.0;
    P_NaK=2.724*1.16;
    alpha_TGF = 0;
elseif Celltype == 3
    G_Ks = 0.098;
    G_to=0.294;
    knaca=1350.0;
    P_NaK=2.724*1.170;
    alpha_TGF = 0;
elseif Celltype == 4
    G_Ks = 0.392;  % nanoS per picoF;
    G_to = 0.073;  % nanoS per picoF;
    knaca=1130.0;
    P_NaK=2.724*1.190;
    alpha_TGF = 0.3;
else
    error('Error:Incorrect input of parameters')
end
% Parameters for I_K1
G_K1=5.405;
% Parameters for I_Na
G_Na=14.838;
% Parameters for I_bNa
G_bNa=0.00029;
% Parameters for I_NaK
K_mK=1.0;
K_mNa=40.0;

% Parameters for I_CaL
G_CaL=0.00003980;
% Parameters for I_bCa
G_bCa=0.000592;
% Parameters for I_NaCa
%     knaca=1000.0;
KmNa_i=87.5;
KmCa=1.38;
gama=0.35;
ksat=0.1;
% Parameters for I_pCa
G_pCa=0.1238;
K_pCa=0.0005;
% Parameters for I_pK
G_pK=0.0146;

%--------------------------------------------------------------------------
%               PARAMETERS FOR INITIAL CONDITIONS
%--------------------------------------------------------------------------
m = 0.001395526321596;
h = 0.770528409603514;
j = 0.770528409602808;
x_r1=2.018256351660449e-04;
x_r2=0.475040393022992;
x_s=0.003212519450373;
r= 5.77774818559522e-08;
s=1;
d=0;
f=1.0;
f_Ca=1.0;
g=1.0;
sd= 2.963079701256671e-05;
sf= 0.999920540500379;
sf2=0.999554105872166;
sfcass= 1.;
sRR= 0.998236073757245;
sOO= 4.057373969085594e-08;

%--------------------------------------------------------------------------
%               PARAMETERS FOR INITIAL  Contraction
%--------------------------------------------------------------------------
TRPNCa = 1.154e-3;
N0=0.998770;
N1=0.367612e-4;
P0=0.112735e-3;
P1=0.148856e-3;
P2=0.408484e-3;
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
V_f = -54.238187021998030;

%--------------------------------------------------------------------------
%             Parameter for simulation duration
%--------------------------------------------------------------------------
dt = 0.005;
xns =5;
disp(['Cycle number: ',num2str(xns)]);
sdur =1.0;
I_stim = 0;
stimstrength=-52.0;  %  stimstrength
stims1 = 1000;
stims2 = 1000;
tbegin=20;
xnstims1 = stims1+tbegin+sdur;%+tbegin;     ???
xnstims2 =stims1*(xns-1)+tbegin+stims2+sdur;%tbegin;    ????
endtime=xnstims2+stims2;  %duration of the simulation
nswitch = 0;
ncounts1=1;
napdd=0;
apdtime = zeros(10,1);

time=0; %timestep(ms)
step=0;

PD =[];  PI =[];
sst=[];ssv=[];ssa=[];ssb=[];ssc=[];ssd=[];sse=[];ssf=[];ssg=[];APV = [];
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
        I_stim = 2*stimstrength;
        skip=1;
        ncount=0;
    end
    
    
    %needed to conpute currents
    E_K=Rtonf*(log(K_o/K_i));
    E_Na=Rtonf*(log(Na_o/Na_i));
    E_Ks=Rtonf*(log((K_o+p_KNa*Na_o)/(K_i+p_KNa*Na_i)));
    E_Ca=0.5*Rtonf*(log(Ca_o/Ca_i));
    Ak1=0.1/(1.+exp(0.06*(v-E_K-200)));
    Bk1=(3.*exp(0.0002*(v-E_K+100))+exp(0.1*(v-E_K-10)))/(1.+exp(-0.5*(v-E_K)));
    rec_iK1=Ak1/(Ak1+Bk1);
    rec_I_NaK=(1.0/(1.0+0.1245*exp(-0.1*v/Rtonf)+0.0353*exp(-v/Rtonf)));
    rec_I_pK=1.0/(1.0+exp((25-v)/5.98));
    
    %Compute currents
    I_Na=G_Na*(m^3)*h*j*(v-E_Na);
    I_CaL=G_CaL*sd*sf*sf2*sfcass*4*(v-15)*(F/Rtonf)...
        *(0.25*exp(2*(v-15)/Rtonf)*CaSS-Ca_o)/(exp(2*(v-15)/Rtonf)-1.);
    
    I_to=G_to*r*s*(v-E_K);
    I_Kr=G_Kr*(sqrt(K_o/5.4))*x_r1*x_r2*(v-E_K);
    I_Ks=G_Ks*x_s*x_s*(v-E_Ks);
    %I_K1=G_K1*x_K1oo*(v-E_K);
    I_K1=G_K1*rec_iK1*(v-E_K);
    I_NaCa=knaca*(1./(KmNa_i*KmNa_i*KmNa_i+Na_o*Na_o*Na_o))*(1./(KmCa+Ca_o))*(1./(1+ksat*exp((gama-1)*v/Rtonf)))*(exp(gama*v/Rtonf)*Na_i*Na_i*Na_i*Ca_o-exp((gama-1)*v*F/(R*T))*Na_o*Na_o*Na_o*Ca_i*2.5);
    I_NaK=P_NaK*(K_o/(K_o+K_mK))*(Na_i/(Na_i+K_mNa))*rec_I_NaK;
    I_pCa=G_pCa*Ca_i/(K_pCa+Ca_i);
    I_pK=G_pK*rec_I_pK*(v-E_K);
    I_bNa=G_bNa*(v-E_Na);
    I_bCa=G_bCa*(v-E_Ca);
    
    %% update concentrationsk
    
    dLTRPN_Ca = k_plus_ltrpn*Ca_i*(LTRPN_tot-LTRPN_Ca)-k_minus_ltrpn*LTRPN_Ca;
    LTRPN_Ca = LTRPN_Ca+dt*dLTRPN_Ca;
    dHTRPN_Ca = k_plus_htrpn*Ca_i*(HTRPN_tot-HTRPN_Ca)-k_minus_htrpn*HTRPN_Ca;
    HTRPN_Ca = HTRPN_Ca+dt*dHTRPN_Ca;
    J_trpn = dLTRPN_Ca+dHTRPN_Ca;
    
    kCaSR=maxsr-((maxsr-minsr)/(1+(EC/CaSR)*(EC/CaSR)));
    k1=k1_/kCaSR;k2=k2_*kCaSR;
    dRR=k4*(1-sRR)-k2*CaSS*sRR;
    sRR=sRR+dt*dRR;
    sOO=k1*CaSS*CaSS*sRR/(k3+k1*CaSS*CaSS);
    Irel=Vrel*sOO*(CaSR-CaSS)*(1-alpha_TGF);
    Ileak=Vleak*(CaSR-Ca_i);
    Iup=Vmaxup/(1.+((Kup*Kup)/(Ca_i*Ca_i)));
    Ixfer=Vxfer*(CaSS-Ca_i);
    CaCSQN=Bufsr*CaSR/(CaSR+Kbufsr);
    dCaSR=dt*(Iup-Irel-Ileak);
    bjsr=Bufsr-CaCSQN-dCaSR-CaSR+Kbufsr;
    cjsr=Kbufsr*(CaCSQN+dCaSR+CaSR);
    CaSR=(sqrt(bjsr*bjsr+4*cjsr)-bjsr)/2;
    CaSSBuf=Bufss*CaSS/(CaSS+Kbufss);
    dCaSS=dt*(-Ixfer*(Vc/Vss)+Irel*(Vsr/Vss)+(-I_CaL*inversevssF2));
    bcss=Bufss-CaSSBuf-dCaSS-CaSS+Kbufss;
    ccss=Kbufss*(CaSSBuf+dCaSS+CaSS);
    CaSS=(sqrt(bcss*bcss+4*ccss)-bcss)/2;
    CaBuf=Bufc*Ca_i/(Ca_i+Kbufc);
    dCa_i=dt*((-(I_bCa+I_pCa-2*I_NaCa)*inverseVcF2)-(Iup+J_trpn-Ileak)*(Vsr/Vc)+Ixfer);
    bc=Bufc-CaBuf-dCa_i-Ca_i+Kbufc;
    cc=Kbufc*(CaBuf+dCa_i+Ca_i);
    Ca_i=(sqrt(bc*bc+4*cc)-bc)/2;
    
    Na_i=Na_i-dt*(I_Na+I_bNa+3*I_NaK+3*I_NaCa)*inverseVcF;
    K_i=K_i-dt*(I_stim+I_K1+I_to+I_Kr+I_Ks-2*I_NaK+I_pK)*inverseVcF;
    
    %--------------------------------------------------------------------------
    %           Compute steady state values and time constants
    %--------------------------------------------------------------------------
    %--------------I_Na--------------%
    alpha_m=1.0/(1.0+exp((-60.0-v)/5.0));
    bate_m=0.1/(1+exp((v+35.0)/5.0))+0.1/(1.0+exp((v-50.0)/200.0));
    t_m=alpha_m*bate_m;
    m_oo=1.0/((1.0+exp((-56.86-v)/9.03))^2);
    if v>=-40.0
        alpha_h=0.0;
        bate_h=0.77/(0.13*(1.0+exp(-(v+10.66)/11.1)));
        alpha_j=0.0;
        bate_j=0.6*exp(0.057*v)/(1.0+exp(-0.1*(v+32.0)));
    else
        alpha_h=0.057*exp(-(v+80.0)/6.8);
        bate_h=2.7*exp(0.079*v)+3.1e5*exp(0.3485*v);
        alpha_j=(-2.5428e4*exp(0.2444*v)-(6.948e-6)*exp(-0.04391*v))*...
            (v+37.78)/(1+exp(0.311*(v+79.23)));
        bate_j=0.02424*exp(-0.01052*v)/(1.0+exp(-0.1378*(v+40.14)));
    end
    t_h=1.0/(alpha_h+bate_h);
    t_j=1.0/(alpha_j+bate_j);
    h_oo=1.0/((1.0+exp((v+71.55)/7.43))^2);
    j_oo=h_oo;
    
    %--------------I_Kr--------------%
    x_r1oo=1.0/(1.0+exp((-26.0-v)/7.0));
    alpha_xr1=450.0/(1.0+exp((-45.0-v)/10.0));
    bate_xr1=6.0/(1.0+exp((v+30.0)/11.5));
    t_xr1=alpha_xr1*bate_xr1;
    
    x_r2oo=1.0/(1.0+exp((v+88.0)/24.0));
    alpha_xr2=3.0/(1+exp((-60.0-v)/20.0));
    bate_xr2=1.12/(1+exp((v-60.0)/20.0));
    t_xr2=alpha_xr2*bate_xr2;
    
    %--------------I_Ks--------------%
    x_soo=1.0/(1+exp((-5.0-v)/14.0));
    alpha_xs=(1400./(sqrt(1.+exp((5.-v)/6))));
    bate_xs=(1./(1.+exp((v-35.)/15.)));
    t_xs=alpha_xs*bate_xs+80;
    
    if Celltype == 1
        r_oo=1./(1.+exp((20-v)/6.));
        s_oo=1./(1.+exp((v+20)/5.));
        t_r=9.5*exp(-(v+40.)*(v+40.)/1800.)+0.8;
        t_s=85.*exp(-(v+45.)*(v+45.)/320.)+5./(1.+exp((v-20.)/5.))+3.;
    elseif Celltype == 2|Celltype == 4
        r_oo=1./(1.+exp((20-v)/6.));
        s_oo=1./(1.+exp((v+28)/5.));
        t_r=9.5*exp(-(v+40.)*(v+40.)/1800.)+0.8;
        t_s=1000.*exp(-(v+67)*(v+67)/1000.)+8.;
    elseif Celltype == 3
        r_oo=1./(1.+exp((20-v)/6.));
        s_oo=1./(1.+exp((v+20)/5.));
        t_r=9.5*exp(-(v+40.)*(v+40.)/1800.)+0.8;
        t_s=85.*exp(-(v+45.)*(v+45.)/320.)+5./(1.+exp((v-20.)/5.))+3.;
    else
        error('Error:Incorrect input of parameters')
    end
    %--------------I_CaL--------------%
    D_INF=1./(1.+exp((-8-v)/7.5));
    Ad=1.4/(1.+exp((-35-v)/13))+0.25;
    Bd=1.4/(1.+exp((v+5)/5));
    Cd=1./(1.+exp((50-v)/20));
    TAU_D=Ad*Bd+Cd;
    F_INF=1./(1.+exp((v+20)/7));
    Af=1102.5*exp(-(v+27)*(v+27)/225);
    Bf=200./(1+exp((13-v)/10.));
    Cf=(180./(1+exp((v+30)/10)))+20;
    TAU_F=Af+Bf+Cf;
    F2_INF=0.67/(1.+exp((v+35)/7))+0.33;
    Af2=600*exp(-(v+25)*(v+25)/170);
    Bf2=31/(1.+exp((25-v)/10));
    Cf2=16/(1.+exp((v+30)/10));
    TAU_F2=Af2+Bf2+Cf2;
    FCaSS_INF=0.6/(1+(CaSS/0.05)*(CaSS/0.05))+0.4;
    TAU_FCaSS=80./(1+(CaSS/0.05)*(CaSS/0.05))+2.;
    
    %Update gates
    m=m_oo-(m_oo-m)*exp(-dt/t_m);
    h=h_oo-(h_oo-h)*exp(-dt/t_h);
    j=j_oo-(j_oo-j)*exp(-dt/t_j);
    x_r1=x_r1oo-(x_r1oo-x_r1)*exp(-dt/t_xr1);
    x_r2=x_r2oo-(x_r2oo-x_r2)*exp(-dt/t_xr2);
    x_s=x_soo-(x_soo-x_s)*exp(-dt/t_xs);
    s=s_oo-(s_oo-s)*exp(-dt/t_s);
    r=r_oo-(r_oo-r)*exp(-dt/t_r);
    sd = D_INF-(D_INF-sd)*exp(-dt/TAU_D);
    sf =F_INF-(F_INF-sf)*exp(-dt/TAU_F);
    sf2 =F2_INF-(F2_INF-sf2)*exp(-dt/TAU_F2);
    sfcass =FCaSS_INF-(FCaSS_INF-sfcass)*exp(-dt/TAU_FCaSS);
    
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
    dN1 = (k_PN*P1-(k_NP+g_10SL)*N1);
    N1=N1+dt*dN1;
    P0=P0+dt*(g_10SL*P1+k_NP*N0-(k_PN+f_01)*P0);
    dP1 = (g_21SL*P2+f_01*P0+k_NP*N1-(f_12+g_10SL+k_PN)*P1);
    P1=P1+dt*dP1;
    dP2 = (g_32SL*P3+f_12*P1-(f_23+g_21SL)*P2);
    P2=P2+dt*dP2;
    dP3= (f_23*P2-g_32SL*P3);
    P3=P3+dt*dP3;
    
    Sigma=g_minxb*2*g_minxb*3*g_minxb+f_01*2*g_minxb*3*g_minxb+f_01*f_12*3*g_minxb+f_01*f_12*f_23;
    P1_max=(f_01*(2*g_minxb)*(3*g_minxb))/(Sigma);
    P2_max=(f_01*f_12*3*g_minxb)/(Sigma);
    P3_max=(f_01*f_12*f_23)/(Sigma);
    F_max=P1_max+2*P2_max+3*P3_max;
    
    F_contrn=-(P1+N1+2*P2+3*P3)/(F_max);
    SL=0.7*F_contrn+SL_0;
    dSL=-0.8*(dP1+dN1+2*dP2+3*dP3)/(F_max);
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
    I_inter_fibro=(V_f-v)*G_gap;
    I_inter_myo= -Num_Fb *I_inter_fibro;
    fI_tot = fI_Na+fI_K1+fI_to+fI_b+fI_NaK+ fI_MGC;
    I_tot=I_Kr+I_Ks+I_K1+I_to+I_Na+I_bNa+I_CaL+I_bCa+I_NaCa+I_NaK+I_pCa+...
        I_pK+I_stim;
    if Num_Fb == 0
        V_f=V_f-dt*fI_tot;
    else
        V_f=V_f-dt*(fI_tot+I_inter_fibro/fC_m);
    end
    
    svol=v-dt*(I_tot+Num_Fb*I_inter_myo/185);%%%%Cm=72;
    
    if(time > (stims1*(xns-2)))&& (time < (stims1*(xns-1)+tbegin))
        APV = [APV; v];
    end
    
    if(time > (stims1*(xns-1)+tbegin))&& (time < (stims1*(xns-1)+sdur+tbegin))
        AP90=max(APV)-(max(APV)-min(APV))*0.9;
    end
    
    if(time > (stims1*(xns-1)+sdur+tbegin))
        if ((svol-AP90)*(v-AP90)<0)
            napdd=napdd+1;
            apdtime(napdd) = time;
        end
        if(napdd==3)
            DI=apdtime(2)-apdtime(1);
            APD=apdtime(3)-apdtime(2);
        end
    end
    v = svol;
    if mod(step,50)==0 % To improve the speed of computer
        sst = [sst;time];
        ssv = [ssv; v];
        ssa = [ssa; V_f];
        ssb = [ssb; F_contr];
        ssc = [ssc; SL];
        ssd = [ssd; Ca_i*1000]; % To change the units
        sse = [sse; Na_i];
        ssf = [ssf; K_i];
    end
    step=step+1;
    time = time+dt;
    
end

AA = [apdtime(1) AP90; apdtime(2) AP90];
BB = [apdtime(2) AP90; apdtime(3) AP90];
if(napdd==3)
    PD =[PD;APD];
    PI =[PI; DI];
end
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
dlmwrite('TNNP_Myofibroblast.txt',M2,'delimiter','\t','precision',6)

t2 = clock;
disp(['total time: ',num2str(etime(t2,t1))]);