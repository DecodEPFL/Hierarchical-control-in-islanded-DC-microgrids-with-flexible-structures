%%% In this script both the network topology and physical parameters are
%%% defined.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Network Data

% n_nodes - number of nodes
% ni - vector of initial nodes of the interconnections
% nf - vector of final nodes of the interconnections

netdata.n_nodes  =                                                     13;

netdata.ni       =    [1,  2,   2,  3,  3,  3,  4,  4,  8,   7,   7,   7];

netdata.nf       =    [2,  3,  11,  5,  6,  4,  7,  8,  9,  10,  12,  13];


% Slack Node

netdata.ns=1;         % Slack Node Position 
netdata.Vs=400;       % Slack Node Nominal Voltage


% Nodes definition

netdata.nodi_gen              =     [5,9,10,11,12,13];                     % Generation Nodes

netdata.nodi_ctrl_gen         =           [5,9,10,11];                     % Fully Controllable Generators' Nodes

netdata.nodi_non_ctrl_gen     =               [12,13];                     % Non Controllable Generators' Nodes

netdata.nodi_load             =          [2,5,6,9,10];                     % Load Nodes

netdata.nodi_non_ctrl         =    [2,5,6,9,10,12,13];                     % Non Controllable Nodes

netdata.gen_batt              =              [5,9,10];                     % Batteries' Nodes
    

% Nominal Frequency and Voltages
netdata.wn=2*pi*50;
netdata.Vn=400;


% Network Parameters defined for each interconnection between ni and nf

% R  - Resistance  [Ohm]
% L  - Inductance  [H]
% C  - Capacitance   [F]
% G  - Conductance [Siemens]
% T  - Elastance   [1/H]
% k  - Trasformer Turns Ratio [V/V]                                                   
        


netdata.R   =    [0.007020;  0.051570;  0.015280;  0.010205;  0.015665;  0.019100;  0.181125;  0.008190;  0.018140; 0.00018140; 0.00018140; 0.00018140];

netdata.L   =    [0.006750;  0.020250;  0.006000;  0.004843;  0.004953;  0.007500;  0.029213;  0.007875;  0.016260; 0.00018140; 0.00016260; 0.00018140]/(2*pi*50);

netdata.C   =    zeros(size(netdata.R));

netdata.G   =    zeros(size(netdata.R));

netdata.T   =    zeros(size(netdata.R));

netdata.k   =    ones(size(netdata.R));                                    % They are all set to 1 since there are no trasformers in the considered microgrid

