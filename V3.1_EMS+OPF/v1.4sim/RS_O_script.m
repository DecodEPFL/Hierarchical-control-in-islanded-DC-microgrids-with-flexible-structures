% %%
% clear all
% clc
% close all
% 

tPLL=2;
t_start_MPC = 10;

T=250;

%% dati carichi

t1=50;
t2=110;
t3=170;

Pbt=1;
Qbt=0;

P1=750e3;
Q1=350e3;

P1v=500e3;
Q1v=250e3;

P2=750e3;
Q2=300e3;

P2v=250e3;
Q2v=250e3;

%% Inizializzazione Ottimizzatore

Dati_rete

n_nodes             =   netdata.n_nodes;
nodi_gen            =   netdata.nodi_gen;
nodi_load           =   netdata.nodi_load;
nodi_ctrl_gen       =   netdata.nodi_ctrl_gen;
nodi_non_crtl_gen   =   netdata.nodi_non_crtl_gen;
n_b                 =   2;

n_states            =   2*n_nodes + 3*n_b + 1;

Dati_MPC
%% Carico profili carichi a potenza costante
S_load         =    zeros( Tsim + N, length(nodi_load));   % Aggiungo N steps di profili in modo da avere la predizione anche negli ultimi minuti

t_step_load_1  = floor((t1-t_start_MPC)/(3600*MPCpar.step))+1;
t_step_load_2  = floor((t2-t_start_MPC)/(3600*MPCpar.step))+1;

S_load(:,1)    =     (P1 + 1i*Q1)*ones(Tsim + N,1);
S_load(t_step_load_1:end,1)  = S_load(t_step_load_1:end,1) + (P1v + 1i*Q1v)*ones(length(S_load(t_step_load_1:end,1)),1);

S_load(:,2)    =     (P2 + 1i*Q2)*ones(Tsim + N,1);
S_load(t_step_load_2:end,2)  = S_load(t_step_load_2:end,2) + (P2v + 1i*Q2v)*ones(length(S_load(t_step_load_2:end,2)),1);

%% Carico profiligeneratori

S_gen = zeros( Tsim + N, length(nodi_gen));
t_step_pv_1   = floor((t3-t_start_MPC)/(3600*MPCpar.step))+1;

S_gen(:,1)                  =     (Ppv1)*ones(Tsim + N,1);
S_gen(t_step_pv_1:end,1)    =     (Ppv2)*ones(length(S_gen(t_step_pv_1:end,1)),1);

%% Definition of the test facility data
% Inital situation for generators and loads reference

P_batt_1 = Pbat1/Srif;
Q_batt_1 = Qbat1/Srif;

P_batt_2 = Pbat2/Srif; 
Q_batt_2 = Qbat2/Srif;

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

var_D = Dprofile(:,2:end)-Dprofile(:,1:end-1); % k+1 -k 
var_D = [zeros(size(var_D,1),1),var_D]./Srif; %k - k-1

%% Definition of initial state

SOC = zeros(3,  1);

SOC_pq(:,1) = [0.8;0.8];
SOC_s(:,1) = 0.6;
SOC(:,1) = [SOC_pq(:,1);SOC_s(:,1)];

%%
var_dist = var_D;
% 
% save('var_dist.mat','var_dist');
