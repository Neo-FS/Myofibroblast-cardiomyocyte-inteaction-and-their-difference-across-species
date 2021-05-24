%Main code from which to run simulations, either from the original model
%of Fusheng et al. (2021)
% Link to Article:
% https://XXXX

clc
clear all
%initial conditions for state variables
V_f = -54.238032147101580;
m_f = 0.166551493609004;
h_f = 0.021549092624088;
j_f = 0.021549092624086;
C_0to = 0.002744277836910;
C_1to = 4.100094085704080e-04;
C_2to = 2.297157828581396e-05;
C_3to = 5.720122491618227e-07;
C_4to = 5.341350664715740e-09;
O_to = 2.284911117674759e-08;

%X0 is the vector for initial sconditions for state variables
X0=[V_f m_f h_f j_f C_0to C_1to C_2to C_3to C_4to O_to]';

CL=1000;%pacing cycle length in ms
beats=1000;%number of beats in the simulation

global tStep tArray
global  fI_Na_store fI_K1_store fI_to_store fI_b_store fI_MGC_store fI_NaK_store Istim_store fI_tot_store

tStep = 1;
tArray = zeros(1,1e7);
fI_Na_store = zeros(1,1e7);
fI_K1_store = zeros(1,1e7);
fI_to_store = zeros(1,1e7);
fI_b_store = zeros(1,1e7);
fI_MGC_store = zeros(1,1e7);
fI_NaK_store = zeros(1,1e7);
Istim_store = zeros(1,1e7);
fI_tot_store = zeros(1,1e7);
%% Run simulation
tic
tspan = [0 50e2];
options = odeset('RelTol',1e-5,'MaxStep',1);
[time,X] = ode15s(@dydt_Fusheng,tspan,X0,options,1);
toc

% %% Output variables
% 
% for n=[1:beats]
%     [time X]=ode15s(@model,[0 CL],X0,options,1);
%     X0=X(size(X,1),:);
%     n; %output beat number to the screen to monitor runtime progress
% end
%% Output variables
tArray = tArray(1:tStep);
fI_Na = fI_Na_store(1:tStep);
fI_K1 = fI_K1_store(1:tStep);
fI_to = fI_to_store(1:tStep);
fI_b = fI_b_store(1:tStep);
fI_NaK = fI_NaK_store(1:tStep);
fI_MGC = fI_MGC_store(1:tStep);
fI_tot = fI_tot_store(1:tStep);
Istim = Istim_store(1:tStep);
fI_K = fI_K1+fI_to;
%rename values in the state variables vector
V_f = X(:,1);
m_f = X(:,2);
h_f = X(:,3);
j_f = X(:,4);
C_0to = X(:,5);
C_1to = X(:,6);
C_2to = X(:,7);
C_3to = X(:,8);
C_4to = X(:,9);
O_to = X(:,10);
%calculate and name dependent variables for the final beat in the
%simulation (i.e. currents and fluxes)

%create plots showing results for the final paced beat
figure (1)
subplot(3,1,1),plot(tArray,Istim),title('fI_tot')
title( 'Simulation protocol,fI_s_t_i_m(pA/pF)','Fontsize',18);
subplot(3,1,2),plot(time,V_f),title('V_f')
title('Membrane potential, V_f(mV)','Fontsize',18);
subplot(3,1,3),plot(tArray,fI_tot),title('fI_tot')
title( 'Totlal current,fI_t_o_t(pA/pF)','Fontsize',18)

figure (2)
subplot(3,2,1),plot(tArray,fI_Na),title('fI_tot')
title( 'Na current,fI_N_a(pA/pF)','Fontsize',18)
subplot(3,2,2),plot(tArray,fI_K1),title('fI_tot')
title( 'Inward rectifying K+ current,fI_K_1(pA/pF)','Fontsize',18)
subplot(3,2,3),plot(tArray,fI_to),title('fI_tot')
title( 'Transient outward K+ current current,fI_t_o(pA/pF)','Fontsize',18)
subplot(3,2,4),plot(tArray,fI_b),title('fI_tot')
title( 'Backup current current,fI_b(pA/pF)','Fontsize',18)
subplot(3,2,5),plot(tArray,fI_MGC),title('fI_tot')
title( 'Mechano-gated current  current,fI_M_G_C(pA/pF)','Fontsize',18)
subplot(3,2,6),plot(tArray,fI_NaK),title('fI_tot')
title( 'NaK exchange current,fI_N_a_K(pA/pF)','Fontsize',18)