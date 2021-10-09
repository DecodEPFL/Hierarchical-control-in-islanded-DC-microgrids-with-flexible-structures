 clear all

load Profiles

%% Network definition

MPCpar                  =   struct;
OPFpar                  =   struct;

load('IncidenceMatrix_DCSlow.mat')

nn = size(B,1);
mm = size(B,2);              % number of lines

sourcenodes         =  1:5;     % Dispatchable generators and batteries
pvnodes             =  6;       % PV systems
loadnodes           =  7:nn;    % ZIP loads
switchnodes         =  [1,4];   % Dispatchable generators

batt_nodes              =   [2,3,5];
MPCpar.nb               =   length(batt_nodes);
nb = MPCpar.nb  ;
MPCpar.batt_nodes      =   batt_nodes;
 
gen_nodes               =   switchnodes;
MPCpar.ng               =   length(gen_nodes);
ng = MPCpar.ng;
MPCpar.gen_nodes       =   gen_nodes;

n  = length(sourcenodes);    % number of generators
p  = length(pvnodes);        % number of PV systems
m  = length(loadnodes);      % number of loads
sn = length(switchnodes);    % number of switchable nodes

% I assume there are three type of loads
realProfiles = struct;

realProfiles.PV(1).YL    =  Profiles.L(1).YL;
realProfiles.PV(1).I     =  Profiles.L(1).I;
realProfiles.PV(1).P     =  1.1*Profiles.L(1).P .*(ones(1,1440) + rand(1,1440)*0.01);

realProfiles.L(1).YL     =  0.07*Profiles.L(2).YL.*(ones(1,1440));
realProfiles.L(1).I      =  0.12*Profiles.L(2).I.*(ones(1,1440) + rand(1,1440)*0.15);
realProfiles.L(1).P      =  0.15*Profiles.L(2).P.*(ones(1,1440) + rand(1,1440)*0.15);
realProfiles.L(2).YL     =  0.07*Profiles.L(2).YL.*(ones(1,1440));
realProfiles.L(2).I      =  0.12*Profiles.L(2).I.*(ones(1,1440) + rand(1,1440)*0.15);
realProfiles.L(2).P      =  0.15*Profiles.L(2).P.*(ones(1,1440) + rand(1,1440)*0.15);
realProfiles.L(3).YL     =  0.07*Profiles.L(2).YL.*(ones(1,1440));
realProfiles.L(3).I      =  0.12*Profiles.L(2).I.*(ones(1,1440) + rand(1,1440)*0.15);
realProfiles.L(3).P      =  0.15*Profiles.L(2).P.*(ones(1,1440) + rand(1,1440)*0.15);

realProfiles.L(4).YL     =  0.06*Profiles.L(3).YL.*(ones(1,1440));
realProfiles.L(4).I      =  0.08*Profiles.L(3).I.*(ones(1,1440) + rand(1,1440)*0.15);
realProfiles.L(4).P      =  0.3*Profiles.L(3).P.*(ones(1,1440) + rand(1,1440)*0.15);
realProfiles.L(5).YL     =  0.06*Profiles.L(3).YL.*(ones(1,1440));
realProfiles.L(5).I      =  0.08*Profiles.L(3).I.*(ones(1,1440) + rand(1,1440)*0.15);
realProfiles.L(5).P      =  0.3*Profiles.L(3).P.*(ones(1,1440) + rand(1,1440)*0.15);
realProfiles.L(6).YL     =  0.06*Profiles.L(3).YL.*(ones(1,1440));
realProfiles.L(6).I      =  0.08*Profiles.L(3).I.*(ones(1,1440) + rand(1,1440)*0.15);
realProfiles.L(6).P      =  0.3*Profiles.L(3).P.*(ones(1,1440) + rand(1,1440)*0.15);

realProfiles.L(7).YL     =  0.08*Profiles.L(4).YL.*(ones(1,1440));
realProfiles.L(7).I      =  0.2*Profiles.L(4).I.*(ones(1,1440) + rand(1,1440)*0.15);
realProfiles.L(7).P      =  0.6*Profiles.L(4).P.*(ones(1,1440) + rand(1,1440)*0.15);
realProfiles.L(8).YL     =  0.08*Profiles.L(4).YL.*(ones(1,1440));
realProfiles.L(8).I      =  0.2*Profiles.L(4).I.*(ones(1,1440) + rand(1,1440)*0.15);
realProfiles.L(8).P      =  0.6*Profiles.L(4).P.*(ones(1,1440) + rand(1,1440)*0.15);
realProfiles.L(9).YL     =  0.08*Profiles.L(4).YL.*(ones(1,1440) );
realProfiles.L(9).I      =  0.2*Profiles.L(4).I.*(ones(1,1440) + rand(1,1440)*0.15);
realProfiles.L(9).P      =  0.6*Profiles.L(4).P.*(ones(1,1440) + rand(1,1440)*0.15);
realProfiles.L(10).YL    =  0.08*Profiles.L(4).YL.*(ones(1,1440) );
realProfiles.L(10).I     =  0.2*Profiles.L(4).I.*(ones(1,1440) + rand(1,1440)*0.15);
realProfiles.L(10).P     =  0.6*Profiles.L(4).P.*(ones(1,1440) + rand(1,1440)*0.15);

delta_nodes = ones(nn,1);     % expresses if a node is connected or not 
delta_lines = ones(mm,1);     % expresses if a line is connected or not 

% ASSUMPTION:  If a generator connected to a line is turned off, the line is disconnected

% Parameters of the filter
C_t     = 2.2*diag(ones(n+p,1))*1e-3;
R_t     = 0.2*diag(ones(n+p,1));
L_t     = 1.9*diag(ones(n+p,1))*1e-3;

% C_t_pv     = 2.2*diag(ones(p,1))*1e-3;
% R_t_pv     = 0.2*diag(ones(p,1));
% L_t_pv     = 1.9*diag(ones(p,1))*1e-3;
% 
% C_l_pv  = 3*diag(ones(p,1))*1e-3;

C_loads = 3*diag(ones(m,1))*1e-3;

%Parameters of the Network
R       = 5*diag(ones(mm,1))*1e-1;
L       = 2*diag(ones(mm,1))*1e-6; 

% For initalization the first instant of the profiles is selected
t = 1;

yl = zeros(1,p+m);
il = zeros(1,p+m);
pl = zeros(1,p+m);

for ii= 1:p
    yl(ii) = [realProfiles.PV(ii).YL(t)];
    il(ii) = [realProfiles.PV(ii).I(t)];
    pl(ii) = [realProfiles.PV(ii).P(t)];
end

for jj= 1:m
    yl(jj+ii) = [realProfiles.L(jj).YL(t)];
    il(jj+ii) = [realProfiles.L(jj).I(t)];
    pl(jj+ii) = [realProfiles.L(jj).P(t)];
end

Y_L     = diag(yl);
I_L     = diag(il);
P_L     = diag(pl);

%Controller Gains 
k1=-eye(n+p);
k2=R_t/30000;
k3=L_t^(-1)*R_t/20;

% new filter parameters 
Ct = blkdiag(C_t,  C_loads);

v_nom  = 100*ones(nn,1);
P_nom = zeros(m+p, length( realProfiles.L(1).I(:)));

for ii= 1:p
    P_nom(ii,:) =  realProfiles.PV(ii).YL(:)'*v_nom(1)^2 + realProfiles.PV(ii).I(:)'*v_nom(1) + realProfiles.PV(ii).P(:)';    
end

for jj = 1:m
    P_nom(p+jj,:) = realProfiles.L(jj).YL(:)'*v_nom(1)^2 + realProfiles.L(jj).I(:)'*v_nom(1) + realProfiles.L(jj).P(:)';     
end

P_nom = [P_nom, P_nom]; % I extend the horizon to avoid problems at the end of the day

PL_nom_s = zeros(m,size(P_nom,2)/15);


for ii = 1 : p
    PL_nom_s(ii,17:81)   =  0.97*smooth(17:81,P_nom(ii,17*15:15:81*15),70,'loess');
    PL_nom_s(ii,:)   =  min( PL_nom_s(ii,:), 0);
end


for jj = 1: m
    PL_nom_s(p+jj,:) = 1*smooth(1:192,P_nom(p+jj,1:15:end),70,'loess');
end
% 
figure
plot(1:15:1440, -PL_nom_s(1,1:96));
hold on
plot(1:1440, -P_nom(1,1:1440));
% % 
figure
plot(1:15:1440, PL_nom_s(5,1:96));
hold on
plot(1:1440, P_nom(5,1:1440));

%% Simulation initialization parameters
%%%%% Initially we assume the PV systems to be pure loads

t_init  = 1;          % seconds for initialization time
dt      = 1e-2;

%%%%%%%%%%%%%%%%%%%%%%%% Initialization
tr = [0:dt:t_init];

p_G_tot      =  sum(P_L*ones(m+p,1)+ diag(v_nom([pvnodes,loadnodes]))*(Y_L)*v_nom([pvnodes,loadnodes]) + (I_L)*v_nom([pvnodes,loadnodes]));        
P_ref_MPC    =  zeros(length(sourcenodes),1);

P_ref_MPC(sourcenodes)     =  p_G_tot/n*ones(1,n);

Delta_P_slack_max      =   ones(n,1)*5000;
Delta_P_slack_min      =  -ones(n,1)*5000;

net_par = struct;
delta_g_k  = ones(1,sn); % I assume all generators ON the beginning
Pred_PV_k = 0;

v_opf = v_nom;

disp('***OPF init started***');
    OPF;
disp('***OPF init finished***');

if  (Pred_PV_k) < 0        
    Vref= [v_opf([sourcenodes pvnodes])];
else
    Vref= [v_opf(sourcenodes);0];
end

t = 0;


% 
% disp('***Simulation init started***');
% tic
% networksimulation
% toc
% disp('***Simulation init finished***');
% 
% figure
% plot(time,yout(:,1:nn))

%% UNITS


rng('default');


MPCpar.SOC_nom          =   0.5;

MPCpar.n_ch             =   [0.9; 0.9; 0.9];
MPCpar.n_dh             =   [0.9; 0.9; 0.9];

MPCpar.Pmax             =   [  80,  50,    40,   80,   60 ];
MPCpar.Pmin             =   [  10, -50,   -40,   10,  -60 ];

MPCpar.Cb               =   [150, 150, 250];                %  batteries capacity

MPCpar.Socmax           =   [ 0.9, 0.9, 0.9 ];
MPCpar.Socmin           =   [ 0.1, 0.1, 0.1 ];

Tsim                    =   24*60;            % Total simulation time in minutes

t_start                 =   1;                  % Starting minute of simulation

dt_sim                  =   60;                 % seconds of simulation iteration for ode 45
dt                      =   0.001;                  % time resolution of simulation in seconds  for ode 45
tr                      =   [0:dt:dt_sim-dt];   % simulation time range for ode45


MPCpar.t_s              =   15;        % MPC sampling time in minutes
OPFpar.t_s              =    3;        % OPF sampling time in minutes


MPCpar.N                =  5*60/MPCpar.t_s;      % MPC prediction horizon

N = MPCpar.N;

MPC_T_iters             =   floor(Tsim/MPCpar.t_s); 


MPCpar.c_b = kron(ones(N,1),    sqrt( 1e-1*ones(nb,1).*( 1 + 0.01*rand(nb,1))));
MPCpar.c_g = kron(ones(N,1),    sqrt( 2*1e2*ones(ng,1).*( 1 + 0.3*rand(ng,1))));

% 
% 
% t_sim                   =   60;        % seconds per simulation
% T                       =   11;        % number of simulation iteration
% Tot_sim                 =   t_sim*T;   % minutes of total simulation time
% 


%% Data saving

% 
% y_plot  = zeros(size(yout,2),Tsim*(length(tr)));
% t_plot  = zeros(size(time,2),Tsim*(length(tr)));
% 

datalog = struct;

datalog.MPC.P_ref     = zeros(n+p,     Tsim);
datalog.MPC.Pred_PV   = zeros(p,       Tsim);
datalog.MPC.sol       = zeros(1,       Tsim);
datalog.MPC.DeltaPmax = zeros(n+p,Tsim);
datalog.MPC.DeltaPmin = zeros(n+p,Tsim);


datalog.OPF.v_opf   = zeros(nn,        Tsim);   
datalog.OPF.Vref    = zeros(n+p,       Tsim);

datalog.OPF.v_opf(:,t_start)   =   v_opf;
datalog.OPF.Vref (:,t_start)   =   Vref;

datalog.OPF.DeltaP    = zeros(n + p,        Tsim);
datalog.OPF.DeltaPmax = zeros(n + p,        Tsim);
datalog.OPF.DeltaPmin = zeros(n + p,        Tsim);

datalog.OPF.Pg       = zeros(n + p,        Tsim);
 
datalog.SOC         = zeros(MPCpar.nb, Tsim);

datalog.OPFsim.v_sim    = zeros(nn,       Tsim);   
datalog.OPFsim.p_G_sim    = zeros(n+p,       Tsim);   


