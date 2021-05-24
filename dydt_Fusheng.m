function output=dydt_Sachse(t,X,flag_ode)
%Main code from which to run simulations, either from the original model
%of Fusheng et al. (2021)

%%  Define model parameters
%physical constants
R=8314.5;  % millijoule_per_mole_kelvin (in membrane) - Ideal gas constant
T=295;     % kelvin (in membrane) - Absolute temperature
F=96487;   % coulomb_per_mole (in membrane) - Faraday constant
Rtonf=(R*T)/F;
C_fm = 18;  % picoF (in membrane) - Membrane capacitance
% extracellular ionic concentrations
fK_o = 5.4;  % millimolar (in standard_ionic_concentrations)
fNa_o = 140.0;  % millimolar (in standard_ionic_concentrations)
% initital intracellular ionic concentrations
fK_i = 129.4349; % millimolar (in standard_ionic_concentrations)
fNa_i = 10;  % millimolar (in standard_ionic_concentrations)

%%   Constants of fibroblast Parameter for [K_o]=5.4mM
GNa_f = 4.2;

P_to = 3.540000000000000e-04;
kvo=30e-3;%
k_vo=2e-3;%
zv=1.28;%
z_v=-1.53;%
ko=77e-3;%
k_o=18e-3;%

fG_K1=0.178;
a_K1=0.94;
b_K1=1.26;

fG_MGC = 0.018;%pS/pF
detaCL = 0;

fG_bNa=4.0e-3;
fG_bK=4.96e-3;

%  NaK
fV_rev=-150.0;
fB=-200;
fK_mk=1.0;%mM3.0
fK_mNa=11.0;%mM
fI_NaKoo=1.75;

%% Global Variable for Time
global tStep tArray
if t > tArray(tStep) % Roughly eliminates data from rejected time steps
    tStep = tStep + 1;
end
tArray(tStep) = t;

% give names to the state vector values
V_f = X(1);
m_f = X(2);
h_f = X(3);
j_f = X(4);
C_0to = X(5);
C_1to = X(6);
C_2to = X(7);
C_3to = X(8);
C_4to = X(9);
O_to = X(10);


%  reversal potentials
fE_Na = Rtonf*log(fNa_o/fNa_i);
fE_K = Rtonf*log(fK_o/fK_i);

%  fI_Na
fm_inf = 1.0/(1.0+exp(-(V_f+42)/7.6));
fh_inf=1.0/(1.0+exp((V_f+84)/7.8));
tm=1.36/((0.32*(V_f+47))/(1.0-exp(-0.1*(V_f+47)))+0.08*exp(-V_f/11));
th = 1.0/(1.03e-5*exp(-(V_f+4.5)/7.5)+0.623*exp((V_f+0.5)/73.5));
tj =1.0/(0.00091*exp(-(V_f+61)/9.6)+0.057*exp((V_f+108)/334));
dm_f = (fm_inf-m_f)/tm;
dh_f = (fh_inf-h_f)/th;
dj_f = (fh_inf-j_f)/tj;
fI_Na = GNa_f*m_f*h_f*j_f*(V_f-fE_Na);
global fI_Na_store
fI_Na_store(tStep) = fI_Na;

%  fI_K1
O_K1=1/(a_K1+exp(b_K1*(V_f-fE_K)/Rtonf));
fI_K1=fG_K1*O_K1*sqrt(fK_o*0.001)*(V_f-fE_K);
global fI_K1_store
fI_K1_store(tStep) = fI_K1;

%  fI_to
kv=kvo*exp(V_f*zv/Rtonf);
k_v=k_vo*exp(V_f*z_v/Rtonf);
if C_0to<0
    C_0to=0;
    dC_0to = 0;
elseif C_0to>1
    C_0to=1;
    dC_0to = 1;
else
    dC_0to=(k_v*C_1to-4*kv*C_0to);
end
if C_1to<0
    C_1to=0;
    dC_1to =0;
elseif C_1to>1
    C_1to=1;
    dC_1to = 1;
else
    dC_1to=(2*k_v*C_2to+4*kv*C_0to-(3*kv+k_v)*C_1to);
end
if C_2to<0
    C_2to=0;
    dC_2to = 0;
elseif C_2to>1
    C_2to=1;
    dC_2to = 0;
else
    dC_2to=(3*k_v*C_3to+3*kv*C_1to-(2*kv+2*k_v)*C_2to);
end
if C_3to<=0
    C_3to=0;
    dC_3to = 0;
elseif C_3to>=1
    C_3to=1;
    dC_3to = 0;
else
    dC_3to=(4*k_v*C_4to+2*kv*C_2to-(kv+3*k_v)*C_3to);
end
if C_4to<=0
    C_4to=0;
    dC_4to = 0;
elseif C_4to>=1
    C_4to=1;
    dC_4to = 0;
else
    dC_4to=(k_o*O_to+kv*C_3to-(ko+4*k_v)*C_4to);
end
if O_to<0
    O_to=0;
    dO_to=0;
elseif O_to>1
    O_to=1;
    dO_to =0;
else
    dO_to=(ko*C_4to-k_o*O_to);
end
fI_to=P_to*O_to*(V_f*F/Rtonf)*(fK_i-fK_o*exp(-V_f/Rtonf))/...
    (1-exp(-V_f/Rtonf));
global fI_to_store
fI_to_store(tStep) = fI_to;

%   fI_MGC
if detaCL <= 0
    fO_MGC = 0;
elseif detaCL >=4.02
    fO_MGC = 1;
else
    fO_MGC = 0.0415*(detaCL^2.286);
end
fI_MGC = fG_MGC*fO_MGC*V_f;
global fI_MGC_store
fI_MGC_store(tStep) = fI_MGC;

%   fI_b
fI_b=fG_bK*(V_f-fE_K)+fG_bNa*(V_f-fE_Na);
global fI_b_store
fI_b_store(tStep) = fI_b;

%   fI_NaK
fI_NaK=fI_NaKoo*(fK_o/(fK_o+fK_mk))*((V_f-fV_rev)/(V_f-fB))*...
    ((fNa_i/(fNa_i+fK_mNa))^1.5);
global fI_NaK_store
fI_NaK_store(tStep) = fI_NaK;

%  simulation time
cycleLength = 500;
stimstrength = -20;
if mod(t-20,cycleLength) <= 1
    Istim = stimstrength;
else
    Istim = 0.0;
end
global Istim_store
Istim_store(tStep) = Istim;

fI_tot = fI_Na+fI_K1+fI_to+fI_b+fI_NaK+ fI_MGC;
global fI_tot_store
fI_tot_store(tStep) = fI_tot;

dV_f = -fI_tot-Istim;

%output the state venctor when ode_flag==1, and the calculated currents and fluxes otherwise
if flag_ode==1
    output=[dV_f dm_f dh_f dj_f dC_0to dC_1to dC_2to dC_3to dC_4to dO_to]';
else
    output=[fI_Na fI_K1 fI_to fI_b fI_MGC fI_NaK];
end

