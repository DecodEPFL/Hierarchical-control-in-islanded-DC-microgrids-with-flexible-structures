MPCpar   = struct;



h_islanded  =  T/3600;
starting_hour = 12;
% Parameters of the cost function
MPCpar.step                 =        (1/6)/60;              % h
MPCpar.N                    =        12;
N                           =        MPCpar.N;

Tsim                        =        floor(h_islanded/MPCpar.step);                % steps

MPCpar.Cf.p_v(1:n_nodes)    =        [1, 1, 1, 1, 1, 1]/(N*0.01^2);
MPCpar.Cf.h                 =        1e5;
MPCpar.Cf.dU                =        [50 50 5 5]/(N*0.2^2);
MPCpar.Cf.Ploss             =        [1e1]/(N*0.005^2);
MPCpar.Cf.P                 =        [1  50 100]/(N*(1)^2);    
MPCpar.Cf.Q                 =        [1 100 100]/(N*(1)^2);   
MPCpar.h_islanded           =        h_islanded;
MPCpar.starting_hour        =        starting_hour;


%% Definition of units capabilities
ft_F=10;

Vrif_M=sqrt(2)/sqrt(3)*20e3;

Vdc_M=800;

Cint=18e-3;

VbatM_min=400;
VbatM_max=1.2*VbatM_min;

VbatM=450;

Re_BpM=4.25e-3;
Ce_BpM=1.30e5;
Ri_BpM=1.00e-1;
Ci_BpM=6.00e4;

ns_BpM=40;
np_BpM=21;

Rext_BpM=Re_BpM*ns_BpM/np_BpM;
Cext_BpM=Ce_BpM*np_BpM/ns_BpM;
Rin_BpM=Ri_BpM*ns_BpM/np_BpM;
Cin_BpM=Ci_BpM*np_BpM/ns_BpM;

SOC_M0=(VbatM-VbatM_min)/(VbatM_max-VbatM_min);

C_batM       = Cext_BpM+Cin_BpM;                           % Battery capacity in farad
C_en_batM_J  = 0.5*C_batM*(VbatM_max^2 - VbatM_min^2);     % Battery available energy in J
C_en_batM    = (1/3600)*C_en_batM_J;                       % Battery available energy in Wh
C_en_batM_pu = C_en_batM/Srif;                             % Battery available energy in puh

%parametri regolatore unidirezionale

P_u=0.005;
I_u=0.1;

%dati 1

Vdc_1=800;

C1=1e-9;

Vbat1_min=400;
Vbat1_max=1.2*Vbat1_min;

Vbat1=450;

Re_Bp1=4.25e-3;
Ce_Bp1=1.30e5;
Ri_Bp1=1.00e-1;
Ci_Bp1=6.00e4;

ns_Bp1=40;
np_Bp1=11;

Rext_Bp1=Re_Bp1*ns_Bp1/np_Bp1;
Cext_Bp1=Ce_Bp1*np_Bp1/ns_Bp1;
Rin_Bp1=Ri_Bp1*ns_Bp1/np_Bp1;
Cin_Bp1=Ci_Bp1*np_Bp1/ns_Bp1;

SOC_10=(Vbat1-Vbat1_min)/(Vbat1_max-Vbat1_min);

Pbat1=150e3;
Qbat1=50e3;

Ppv1=500e3;
Ppv2=850e3;

C_bat1       = Cext_Bp1+Cin_Bp1;                           % Battery capacity in farad
C_en_bat1_J  = 0.5*C_bat1*(Vbat1_max^2 - Vbat1_min^2);     % Battery available energy in J
C_en_bat1    = (1/3600)*C_en_bat1_J;                    % Battery available energy in Wh
C_en_bat1_pu = C_en_bat1/Srif;                       % Battery available energy in puh

%dati 2

Vdc_2=800;

C2=1e-9;

Vbat2_min=400;
Vbat2_max=1.2*Vbat2_min;

Vbat2=450;

Re_Bp2=4.25e-3;
Ce_Bp2=1.30e5;
Ri_Bp2=1.00e-1;
Ci_Bp2=6.00e4;

ns_Bp2=40;
np_Bp2=11;

Rext_Bp2=Re_Bp2*ns_Bp2/np_Bp2;
Cext_Bp2=Ce_Bp2*np_Bp2/ns_Bp2;
Rin_Bp2=Ri_Bp2*ns_Bp2/np_Bp2;
Cin_Bp2=Ci_Bp2*np_Bp2/ns_Bp2;

SOC_20=(Vbat2-Vbat2_min)/(Vbat2_max-Vbat2_min);

Pbat2=250e3;
Qbat2=50e3;

C_bat2       = Cext_Bp2+Cin_Bp2;                           % Battery capacity in farad
C_en_bat2_J  = 0.5*C_bat2*(Vbat2_max^2 - Vbat2_min^2);     % Battery available energy in J
C_en_bat2    = (1/3600)*C_en_bat2_J;                    % Battery available energy in Wh
C_en_bat2_pu = C_en_bat2/Srif;                       % Battery available energy in puh

%Parametri PLL%

k=sqrt(2);

Vd=Vrif_M;

ft_PLL=50;

w_t=2*pi*10;

kp_PLL_max=w_t/Vd;
kp_PLL=0.9*kp_PLL_max;
ki_PLL=w_t/Vd*(-w_t+sqrt(2*w_t^2-Vd^2*kp_PLL^2));


% Vincoli slack generator

MPCpar.Pmax_ns_pu   =  2;
MPCpar.Pmin_ns_pu   = -2;
MPCpar.Qmax_ns_pu   =  1;
MPCpar.Qmin_ns_pu   = -1;
MPCpar.Socmax_ns    =  1;
MPCpar.Socmin_ns    =  0;
MPCpar.C_ns         =  C_en_batM_pu;  %puh   

% Vincoli gen 1
MPCpar.Pmax_gen1_pu =  1;
MPCpar.Pmin_gen1_pu = -1;
MPCpar.Qmax_gen1_pu =  1;
MPCpar.Qmin_gen1_pu = -1;
MPCpar.Socmax_gen1  =  1;
MPCpar.Socmin_gen1  =  0;
MPCpar.C_gen1       =  C_en_bat1_pu; %puh

% Vincoli gen 2
MPCpar.Pmax_gen2_pu =  1;
MPCpar.Pmin_gen2_pu = -1;
MPCpar.Qmax_gen2_pu =  1;
MPCpar.Qmin_gen2_pu = -1;
MPCpar.Socmax_gen2  =  1;
MPCpar.Socmin_gen2  =  0;
MPCpar.C_gen2       =  C_en_bat2_pu; %puh

% Vincoli load rete
MPCpar.MaxPLoadStep_up_pu  =  0.5;
MPCpar.MaxPLoadStep_dwn_pu =  0.5;
MPCpar.MaxQLoadStep_up_pu  =  0.5;
MPCpar.MaxQLoadStep_dwn_pu =  0.5;

MPCpar.Vmax         =  1.1;
MPCpar.Vmin         =  0.9;
