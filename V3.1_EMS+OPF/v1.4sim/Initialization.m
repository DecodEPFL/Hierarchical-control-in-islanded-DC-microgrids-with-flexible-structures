% Initialization

Dati_rete

n_nodes             =   netdata.n_nodes;
nodi_gen            =   netdata.nodi_gen;
nodi_load           =   netdata.nodi_load;
nodi_ctrl_gen       =   netdata.nodi_ctrl_gen;
nodi_non_crtl_gen   =   netdata.nodi_non_crtl_gen;
% nodi_non_crtl     =   netdata.nodi_non_crtl;
n_b                 =   2;

n_states            =   2*n_nodes + 3*n_b + 1;

%% Definition of the control parameters
MPCpar   = struct;

h_islanded  =  60/60;
starting_hour = 12;
% Parameters of the cost function
MPCpar.step                 =        (1/3)/60;              % h
MPCpar.N                    =        9;
N=MPCpar.N  ;

Tsim                        =        floor(h_islanded/MPCpar.step);                % steps

MPCpar.Cf.p_v(1:n_nodes)    =        [1, 1, 1, 1, 1, 1]/(N*0.01^2);
MPCpar.Cf.h                 =        1e5;
MPCpar.Cf.dU                =        1e1*[50 50 5 5]/(N*0.2^2);
MPCpar.Cf.Ploss             =        [1e1]/(N*0.005^2);
MPCpar.Cf.P                 =        [1 50 100]/(N*(1)^2);    
MPCpar.Cf.Q                 =        [1 100 100]/(N*(1)^2);    

%% Carico profili carichi a potenza costante
S_load      =       zeros( Tsim + N, length(nodi_load));   % Aggiungo N steps di profili in modo da avere la predizione anche negli ultimi minuti

load Load_profiles
load P_profile

P_load_1 =  P_profile(:,2)*20;
Q_load_1 =  P_profile(:,2)*20*tan(acos(0.8));

S_load(:,1) =  (interp1(1:1:(h_islanded)*60 + N*MPCpar.step*60 ,P_load_1(starting_hour*60:1:(starting_hour+h_islanded)*60+N*MPCpar.step*60 -1),MPCpar.step*60:MPCpar.step*60:(h_islanded)*60 + N*MPCpar.step*60,'linear','extrap')+ 1i*interp1(1:1:(h_islanded)*60+ N*MPCpar.step*60,Q_load_1(starting_hour*60:1:(starting_hour+h_islanded)*60+ N*MPCpar.step*60-1),MPCpar.step*60:MPCpar.step*60:(h_islanded)*60+ N*MPCpar.step*60,'linear','extrap'));
S_load(:,2) =  (interp1(1:1:(h_islanded)*60 + N*MPCpar.step*60 ,P_load_2(starting_hour*60:1:(starting_hour+h_islanded)*60+N*MPCpar.step*60 -1),MPCpar.step*60:MPCpar.step*60:(h_islanded)*60 + N*MPCpar.step*60,'linear','extrap')+ 1i*interp1(1:1:(h_islanded)*60+ N*MPCpar.step*60,Q_load_2(starting_hour*60:1:(starting_hour+h_islanded)*60+ N*MPCpar.step*60-1),MPCpar.step*60:MPCpar.step*60:(h_islanded)*60+ N*MPCpar.step*60,'linear','extrap'));

% S_load(:,1) =  400e3*ones(size(S_load(:,1),1),1) + 1i*100e3*ones(size(S_load(:,1),1),1);
% S_load(:,2) =  200e3*ones(size(S_load(:,1),1),1) + 1i*100e3*ones(size(S_load(:,1),1),1);

S_load(30:50,2) = S_load(30:50,2) + 0.5*Srif*ones(21,1);
S_load(120:150,2) = S_load(120:150,2) + 0.5*Srif*ones(31,1);
S_load(70:95,2) = S_load(70:95,2) + 1i*0.5*Srif*ones(26,1);
%% Carico profiligeneratori

S_gen = zeros( Tsim + N, length(nodi_gen));

load P_solar  % Nodo 5

S_gen(:,1) =  5*(interp1(1:1:(h_islanded)*60 + N*MPCpar.step*60,Psolar(starting_hour*60:1:(starting_hour+h_islanded)*60-1 + N*MPCpar.step*60)',MPCpar.step*60:MPCpar.step*60:(h_islanded)*60 + N*MPCpar.step*60,'linear','extrap'));
% S_gen(:,1) =  100e3*ones(size(S_load(:,1),1),1);

%% Definition of the test facility data
% Inital situation for generators and loads reference

P_batt_1 = 200e3/Srif;
Q_batt_1 = 100e3/Srif;

P_batt_2 = 200e3/Srif; 
Q_batt_2 = 100e3/Srif;

% parametri regolatore generatore PQ 1
RegG(1).nG = 5;
RegG(1).P  = P_batt_1 + real(S_gen(1,1))/Srif;
RegG(1).Q  = Q_batt_1;

% parametri regolatore generatore PQ 2
RegG(2).nG = 6;
RegG(2).P  = P_batt_2;
RegG(2).Q  = Q_batt_2;

% parametri carico C 2
RegL(1).nL =  2;
RegL(1).P  = -real(S_load(1,1))/Srif;
RegL(1).Q  = -imag(S_load(1,1))/Srif;

% parametri carico C 1
RegL(2).nL =  3;
RegL(2).P  = -real(S_load(1,2))/Srif;
RegL(2).Q  = -imag(S_load(1,2))/Srif;

%% Definition of units capabilities

% Vincoli slack generator
MPCpar.Pmax_ns_pu   =  2;
MPCpar.Pmin_ns_pu   = -2;
MPCpar.Qmax_ns_pu   =  1;
MPCpar.Qmin_ns_pu   = -1;
MPCpar.Socmax_ns    =  1;
MPCpar.Socmin_ns    =  0;
MPCpar.C_ns         =  1;  %puh   

% Vincoli gen 1
MPCpar.Pmax_gen1_pu =  1;
MPCpar.Pmin_gen1_pu = -1;
MPCpar.Qmax_gen1_pu =  1;
MPCpar.Qmin_gen1_pu = -1;
MPCpar.Socmax_gen1  =  1;
MPCpar.Socmin_gen1  =  0;
MPCpar.C_gen1       =  0.5; %puh

% Vincoli gen 2
MPCpar.Pmax_gen2_pu =  1;
MPCpar.Pmin_gen2_pu = -1;
MPCpar.Qmax_gen2_pu =  1;
MPCpar.Qmin_gen2_pu = -1;
MPCpar.Socmax_gen2  =  1;
MPCpar.Socmin_gen2  =  0;
MPCpar.C_gen2       =  0.5; %puh


% Vincoli load rete
MPCpar.MaxPLoadStep_up_pu  =  0.5;
MPCpar.MaxPLoadStep_dwn_pu =  0.5;
MPCpar.MaxQLoadStep_up_pu  =  0.1;
MPCpar.MaxQLoadStep_dwn_pu =  0.1;

MPCpar.Vmax         =  1.1;
MPCpar.Vmin         =  0.9;

%% Definition of disturbance profiles

% Definition of disturbance of reference powers for each node; transit nodes
% have no disturbance since there are neither generetors nor loads.

Dprofile      =       zeros((n_nodes-1)*2, Tsim + N); %The slack node is not affected by disturbances

%%% Definition of non controllable generators disturbance (solar panels and wind turbines)

Dprofile (1,:)   =   -real(S_load(:,1));              
Dprofile (1 + n_nodes-1,:)   =   -imag(S_load(:,1));      
 
Dprofile (2,:)   =   -real(S_load(:,2));              
Dprofile (2+ n_nodes-1,:)   =   -imag(S_load(:,2));    

Dprofile (4,:)   =    real(S_gen(:,1));              
Dprofile (4 + n_nodes-1,:)   =    imag(S_gen(:,1));      


%%% Definition of loads disturbance

var_D = Dprofile(:,2:end)-Dprofile(:,1:end-1);
var_D = [var_D,var_D(:,end)]./Srif;

%% Definition of initial state

SOC = zeros(3,  1);

SOC_pq(:,1) = [0.8;0.8];
SOC_s(:,1) = 0.6;
SOC(:,1) = [SOC_pq(:,1);SOC_s(:,1)];

%% DataLog


datalog = struct;

datalog.P_loss   = zeros(n_nodes-1,Tsim);

datalog.P_batt_1 = zeros(1,Tsim);
datalog.P_batt_2 = zeros(1,Tsim);   
datalog.Q_batt_1 = zeros(1,Tsim);
datalog.Q_batt_2 = zeros(1,Tsim);
datalog.P_node_5 = zeros(1,Tsim+N);

datalog.P_batt_1(1) = P_batt_1;
datalog.P_node_5(1) = P_batt_1 + real(S_gen(1,1))/Srif;
datalog.P_batt_2(1) = P_batt_2;   
datalog.Q_batt_1(1) = Q_batt_1;
datalog.Q_batt_2(1) = Q_batt_2;

datalog.Umpc   =  zeros(length(nodi_ctrl_gen)*2,N,Tsim);
datalog.Jopt   =  zeros(1,Tsim);
datalog.h      =  zeros(n_nodes-1+1+1,N,Tsim);
    
datalog.x      =  zeros(n_states,Tsim);
datalog.SOC    =  zeros(3,Tsim);
datalog.SOC(:,1) = SOC(:,1);