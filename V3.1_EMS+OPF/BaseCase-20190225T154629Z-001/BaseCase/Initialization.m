%% Initialization

%% Network data

TestFacility_RealData                                                      % Network Data for the simulator (Real Parameters)
      

n_nodes             =   netdata.n_nodes;
nodi_gen            =   netdata.nodi_gen;
nodi_load           =   netdata.nodi_load;
nodi_ctrl_gen       =   netdata.nodi_ctrl_gen;
nodi_non_ctrl_gen   =   netdata.nodi_non_ctrl_gen;
nodi_non_ctrl       =   netdata.nodi_non_ctrl;
gen_batt            =   netdata.gen_batt;


max_simulation_time = 24;                                                  % Maximum number of simulated hours
n_min = max_simulation_time*60;                                            % Maximum number of simulated minutes

%% Loads Power Profiles - Defined at each MINUTE

Pload=zeros(n_min,length(nodi_load));                                      % Complex Load Power inizialization


%%% Costant Loads's Power Profiles

Cp0  =  [ 4.50;   0.19;   0.28]*1000;                                      % Active Power [W]
Cq0  =  [-4.96;  -2.12;  -0.83]*1000;                                      % Reactive Power [Var]

Pload(:,2)  =    Cp0(1)*ones(n_min,1) + 1i*Cq0(1)*ones(n_min,1);           % Node 5
Pload(:,4)  =    Cp0(2)*ones(n_min,1) + 1i*Cq0(2)*ones(n_min,1);           % Node 7
Pload(:,5)  =    Cp0(3)*ones(n_min,1) + 1i*Cq0(3)*ones(n_min,1);           % Node 9


%%% Varying Loads' Power Profiles 

load Pr_profile                                                            % Active    Resistive   Load  Power

load Qind_profile                                                          % Reactive  Inductive   Load  Power
load Qcap_profile                                                          % Reactive  Capacitive  Load  Power

Pload(:,1) = Pr_profile(:,1) + 1i*(Qind_profile(:,1) + Qcap_profile(:,1)); % Node 2
Pload(:,3) = Pr_profile(:,2) + 1i*(Qind_profile(:,2) + Qcap_profile(:,2)); % Node 6


%% Generators' Parameters

%%% Generators Power Profiles - Defined at each MINUTE

Pgen= zeros(n_min,length(nodi_gen));

%%% Uncontrollable Generators' Power Profiles

load Psolar
load Pwindturbine

%%% Initial Power References for Controllable Generators 

Pgen(:,1)   =    (1000)*ones(n_min,1)     +    1i*(0)*ones(n_min,1);
Pgen(:,2)   =    (1000)*ones(n_min,1)     +    1i*(0)*ones(n_min,1);
Pgen(:,3)   =    (1000)*ones(n_min,1)     +    1i*(0)*ones(n_min,1);
Pgen(:,4)   =    (1000)*ones(n_min,1)     +    1i*(0)*ones(n_min,1);

%%% Power Profiles for Renewable Sources

Pgen(:,5)   =    Psolar;
Pgen(:,6)   =    Pwindturbine;    

%%% Generators Power Limits

% The Active Power Limits are defined for each controllable generator

Gen_Pmax0        =   [  30000,    25000,   12000,   50000];                 
Gen_Pmin0        =   [ -11000,   -25000,   -9000,       0]; 


% The Reactive Power Limits are defined for all generators ( they depends
% on the inverters)

Gen_Qmax0        =   [  60000,    15000,   13000,   38553,    10000,   10000];     
Gen_Qmin0        =   [ -60000,   -15000,  -11000,  -38553,   -10000,  -10000]; 


% The power limits are initially defined through two equal set of variables
% (e.g. Gen_Pmax and Gen_Pmax0)
% However Gen_Pmax might be modified in some cases; for instance
% if a battery is discharged, the maximum output power is then set to zero
% so that the battery can only absorb power.

Gen_Pmax        =   Gen_Pmax0;      Gen_Pmin  =   Gen_Pmin0; 
Gen_Qmax        =   Gen_Qmax0;      Gen_Qmin  =   Gen_Qmin0; 


% Initial value of the frequency integrator output
v               =      0;                                                                 


% Initial states of charge of batteries 

Soc             =     [80;80;80];  


%% Interpolation of Power Profiles at the simulation time step

n=2;                                                                       % Number of defined simulation steps within one minute
n_sim=n*n_min;                                                             % Maximum number of simulation steps

% An additional minute is added to the power profiles in order to perform
% the interpolation also for the last step

Pgen            =      [Pgen;Pgen(end,:)];
Pload           =      [Pload;Pload(end,:)];
Pr_profile      =      [Pr_profile;Pr_profile(end,:)];
Qind_profile    =      [Qind_profile;Qind_profile(end,:)];
Qcap_profile    =      [Qcap_profile;Qcap_profile(end,:)];
Psolar          =      [Psolar;Psolar(end,:)];
Pwindturbine    =      [Pwindturbine;Pwindturbine(end,:)];

% The interpolation process is then performed
Pgen            =      interp1(    1:(n_min+1),    Pgen,           ( 1 : 1/n : (n_min+1))  );
Pload           =      interp1(    1:(n_min+1),    Pload,          ( 1 : 1/n : (n_min+1))  );
Pr_profile      =      interp1(    1:(n_min+1),    Pr_profile,     ( 1 : 1/n : (n_min+1))  );
Qind_profile    =      interp1(    1:(n_min+1),    Qind_profile,   ( 1 : 1/n : (n_min+1))  );
Qcap_profile    =      interp1(    1:(n_min+1),    Qcap_profile,   ( 1 : 1/n : (n_min+1))  );
Psolar          =      interp1(    1:(n_min+1),    Psolar,         ( 1 : 1/n : (n_min+1))  );
Pwindturbine    =      interp1(    1:(n_min+1),    Pwindturbine,   ( 1 : 1/n : (n_min+1))  );

% The last time step is then deleted in order to have exactly n_sim steps
Pgen            =      Pgen(1:end-1,:);
Pload           =      Pload(1:end-1,:);
Pr_profile      =      Pr_profile(1:end-1,:);
Qind_profile    =      Qind_profile(1:end-1,:);
Qcap_profile    =      Qcap_profile(1:end-1,:);
Psolar          =      Psolar(1:end-1,:);
Pwindturbine    =      Pwindturbine(1:end-1,:);

%% Definition of the secondary controller parameters

MPCpar          =   struct;


%%% CONSTRAINTS

    % Bound for voltages
    MPCpar.Vmax     =   440*ones(n_nodes,1);                                    
    MPCpar.Vmin     =   360*ones(n_nodes,1);

    % Bounds for Soc

        % Hard Constraints
    MPCpar.SocMax   =   100*ones(length(gen_batt),1);
    MPCpar.SocMin   =     0*ones(length(gen_batt),1);

        % Soft Constraints
    MPCpar.SocUB    =    80*ones(length(gen_batt),1);
    MPCpar.SocLB    =    20*ones(length(gen_batt),1);

    % Bounds for Control Inputs 
    uvarnom=3000;                                                              % [W]
    MPCpar.uvarnom=uvarnom;

    %%%MPC Horizon
    MPCpar.N=5;
    N=MPCpar.N;

    %Number of times that MPC is performed
    MPCpar.T=n_sim;                                                         
    T  =   MPCpar.T;


%%% COST FUNCTION

% Weights
    MPCpar.Cf.ro                  =                   10000/(10);              % slack variable of voltages
    MPCpar.Cf.rovar               =                  1000/(3000);              % slack variable of reference powers (e.g. in case of battery discharge) 
    MPCpar.Cf.rosoc               =                     3000/(1);              % slack variable for SOCs

    MPCpar.Cf.p_int               =                         1;               % frequency integrator state

    MPCpar.Cf.p_v                 =           10/((1^2)*n_nodes*N);             % voltage difference wrt nominal value

    MPCpar.Cf.p_u1_p              =        0.1/(Gen_Pmax(1)^2*N);              % reference active power of 1st generator
    MPCpar.Cf.p_u1_q              =         10/(Gen_Qmax(1)^2*N);              % reference reactive power of 1st generator

    MPCpar.Cf.p_u2_p              =        0.1/(Gen_Pmax(2)^2*N);              % reference active power of 2nd generator
    MPCpar.Cf.p_u2_q              =         10/(Gen_Qmax(2)^2*N);              % reference reactive power of 2nd generator

    MPCpar.Cf.p_u3_p              =        0.1/(Gen_Pmax(3)^2*N);              % reference active power of 3rd generator
    MPCpar.Cf.p_u3_q              =         10/(Gen_Qmax(3)^2*N);              % reference reactive power of 3rd generator

    MPCpar.Cf.p_u4_p              =        0.1/(Gen_Pmax(4)^2*N);              % reference  active power of 4th generator
    MPCpar.Cf.p_u4_q              =         10/(Gen_Qmax(4)^2*N);              % reference reactive power of 4th generator


    MPCpar.Cf.p_u1_dp             =            100/(uvarnom^2*N);              % reference active power of 1st generator
    MPCpar.Cf.p_u1_dq             =             10/(uvarnom^2*N);              % reference reactive power of 1st generator

    MPCpar.Cf.p_u2_dp             =            100/(uvarnom^2*N);              % reference active power of 2nd generator
    MPCpar.Cf.p_u2_dq             =             10/(uvarnom^2*N);              % reference reactive power of 2nd generator

    MPCpar.Cf.p_u3_dp             =            100/(uvarnom^2*N);              % reference active power of 3rd generator
    MPCpar.Cf.p_u3_dq             =             10/(uvarnom^2*N);              % reference reactive power of 3rd generator

    MPCpar.Cf.p_u4_dp             =            100/(uvarnom^2*N);              % reference  active power of 4th generator
    MPCpar.Cf.p_u4_dq             =             10/(uvarnom^2*N);              % reference reactive power of 4th generator



%% Definition of disturbance profiles

% Definition of disturbance of reference powers for each node; transit nodes
% have no disturbance since there are neither generetors nor loads.

Dprofile      =       zeros(n_nodes*2,T+N);

%%% Definition of non controllable generators disturbance (solar panels and wind turbines)

Dprofile (12,           1:T)   =   real(Pgen(1:T,5)');              
Dprofile (12 + n_nodes, 1:T)   =   imag(Pgen(1:T,5)');      
 
Dprofile (13,           1:T)   =   real(Pgen(1:T,6)');              
Dprofile (13 + n_nodes, 1:T)   =   imag(Pgen(1:T,6)');    

%%% Definition of loads disturbance

for i = 1:length(nodi_load)
    
for j=1:n_nodes
    
    if (j == nodi_load(i))
       
        Dprofile(j,1:n_sim)            =   Dprofile(j,1:n_sim)              -  real(Pload(1:n_sim,i)');
        Dprofile(j+n_nodes, 1:n_sim)   =   Dprofile(j+n_nodes, 1:n_sim)     -  imag(Pload(1:n_sim,i)');
        
    end
    
end

end

Dprofile(:,(n_sim+1):(n_sim+N))=repmat(Dprofile(:,n_sim),1,N);

%%% Definition of disturbance variation so that it is possible to give it to MPC 
% var_D(k) = Dprofile(k)-Dprofile(k-1);

var_D   = zeros(n_nodes*2,n_sim+N);

for i = 1:(n_sim+N-1)
    
var_D(:,i+1) = Dprofile(:,i+1) - Dprofile(:,i);

end

%% Noise is added to the load power profiles

noise_Pr       =    500*randn(T,2);
noise_Qind     =    100*randn(T,2);
noise_Qcap     =    100*randn(T,2);


for i=1:2:T
   
Pr_profile(i,:)       =      Pr_profile(i,:)     +        noise_Pr(i,:);
Qind_profile(i,:)     =      Qind_profile(i,:)   +        noise_Qind(i,:);
Qcap_profile(i,:)     =      Qcap_profile(i,:)   +        noise_Qcap(i,:);

end

% The adding of random noise may cause the absorbed active power to take values
% lower than zero, therefore only positive values are taken

Pr_profile=max(Pr_profile,0);

% Having defined the new power profiles now the Pload vector is again
% defined

Pload(:,1) = Pr_profile(:,1) + 1i*(Qind_profile(:,1) + Qcap_profile(:,1)); % Node 2
Pload(:,3) = Pr_profile(:,2) + 1i*(Qind_profile(:,2) + Qcap_profile(:,2)); % Node 6


%% Definition of initial conditions for the first simulation step

% Now the power entities for each unit are defined for the first simulation
% step, these will be updated at each simulation iteration


Pr1             =   Pr_profile(1,1);                                       % Active Resistive Power of Load 1
Pr2             =   Pr_profile(1,2);                                       % Active Resistive Power of Load 3

Qind1           =   Qind_profile(1,1);                                     % Reactive Inductive Power of Load 1
Qind2           =   Qind_profile(1,2);                                     % Reactive Inductive Power of Load 3

Qcap1           =   Qcap_profile(1,1);                                     % Reactive Inductive Power of Load 1
Qcap2           =   Qcap_profile(1,2);                                     % Reactive Inductive Power of Load 3

SL2             =   Pload(1,2);                                            % Complex Power of Load 2
SL4             =   Pload(1,4);                                            % Complex Power of Load 4
SL5             =   Pload(1,5);                                            % Complex Power of Load 5


% The reference powers are defined for each generation unit. Also these are 
% updated at each simulation step depending on the control action and on the 
% renewable power profiles.

Gen_Pref        =   real([Pgen(1,1), Pgen(1,2), Pgen(1,3), Pgen(1,4), Pgen(1,5), Pgen(1,6)]);
Gen_Qref        =   imag([Pgen(1,1), Pgen(1,2), Pgen(1,3), Pgen(1,4), Pgen(1,5), Pgen(1,6)]);
                                              

% The following script defines the generators' and loads' characteristics for the simulator

generatori_e_carichi;

%% DataLog

% A struct is created so that it is possible to save the simulation data

datalog             =   struct;

datalog.DU          =   zeros(length(nodi_ctrl_gen)*2,T);                  % variation of reference power, MPC control unit

datalog.U           =   zeros(length(nodi_ctrl_gen)*2,T);                  % Total reference powers imposed by MPC

datalog.U(1:4,1)    =   Gen_Pref(1,1:4)';                                  % Reference Active Powers
datalog.U(5:8,1)    =   Gen_Qref(1,1:4)';                                  % Reference Reactive Powers

datalog.X           =   zeros(n_nodes*2+1+length(gen_batt),T);             % State Variables

datalog.Jopt        =   zeros(1,T);                                        % Cost Function

datalog.hsoft       =   zeros(1,T);                                        % slack variable for voltage bounds

datalog.hpsoft      =   zeros(1,T);                                        % slack variable for reference power variations

datalog.Socmpc      =   zeros(length(gen_batt),T);                         % SOC estimation by MPC

datalog.Sgen        =   zeros(length(nodi_gen),T);                         % real generated powers
datalog.Sload       =   zeros(length(nodi_load),T);                        % real absorbed powers

datalog.Gen_Pmax    =   zeros(length(nodi_ctrl_gen),T);                    % maximum active power
datalog.Gen_Pmin    =   zeros(length(nodi_ctrl_gen),T);                    % minimum reactive power 

datalog.Pgen        =   zeros(length(nodi_ctrl_gen)*2,T);                  % real generated power estimation by MPC

