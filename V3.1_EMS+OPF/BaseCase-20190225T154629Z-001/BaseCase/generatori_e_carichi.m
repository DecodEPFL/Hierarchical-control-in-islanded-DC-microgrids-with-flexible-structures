%%% In this script the parameters of all microgrid units are defined

%% GENERATOR 1, N0DE 5, LITHIUM BATTERY

Ctot(1)  = 32000;                                                           % Battery Capacity [Wh]
coeff(1) = 1;                                                               % Charge/Discharge Coefficient

k(1) = 100*coeff(1)*1/(60*n*Ctot(1));                                            

RegG(1).nG = 5;                                                              % Generator Node

% The reference active and reactive power at each time step are defined
RegG(1).Gen_Qref = Gen_Qref(1);                           
RegG(1).Gen_Pref = Gen_Pref(1);

% Definition of Droop Parameters

RegG(1).V = [0.9,0.97,1.03,1.10]*400;         
RegG(1).w = 2*pi*50*[0.98,0.994,1.006,1.02];

RegG(1).Qmax = Gen_Qmax(1);
RegG(1).Qmin = Gen_Qmin(1);

RegG(1).Pmax = Gen_Pmax(1);
RegG(1).Pmin = Gen_Pmin(1);

    % Droop Parameters

RegG(1).mv =  - (  (RegG(1).Pmax)    -   (RegG(1).Pmin)   )   /   (    (   RegG(1).V(4)   -   RegG(1).V(1)  )     );      
RegG(1).mf =    (  (RegG(1).Qmax)    -   (RegG(1).Qmin)   )   /   (    (   RegG(1).w(4)   -   RegG(1).w(1)  )     );

%% GENERATOR 2, NODE 9, LOCCIONI BATTERY

Ctot(2)  = 30000;                                                          % Battery Capacity [Wh]
coeff(2) = 1;                                                              % Charge/Discharge Coefficient

k(2) = 100*coeff(2)*1/(60*n*Ctot(2));        

RegG(2).nG = 9;                                                            % Generator Node

% The reference active and reactive power at each time step are defined
RegG(2).Gen_Qref = Gen_Qref(2);
RegG(2).Gen_Pref = Gen_Pref(2);

% Definition of Droop Parameters
RegG(2).V = [0.9,0.97,1.03,1.10]*400;
RegG(2).w = 2*pi*50*[0.98,0.994,1.006,1.02];

RegG(2).Qmax = Gen_Qmax(2);
RegG(2).Qmin = Gen_Qmin(2);

RegG(2).Pmax = Gen_Pmax(2);
RegG(2).Pmin = Gen_Pmin(2);

    % Droop Parameters
    
RegG(2).mv = -(  (RegG(2).Pmax)    -   (RegG(2).Pmin) )   / ((RegG(2).V(4)-RegG(2).V(1)));        
RegG(2).mf =  (  (RegG(2).Qmax)    -   (RegG(2).Qmin) )   / ((RegG(2).w(4)-RegG(2).w(1)));

%% GENERATOR 3, NODE 10, LEAD BATTERY

Ctot(3) = 55000;                                                           % Battery Capacity [Wh]
coeff(3) = 1;                                                              % Charge/Discharge Coefficient

k(3) = 100*coeff(3)*1/(60*n*Ctot(3));      

RegG(3).nG = 10;                                                           % Generator Node

% The reference active and reactive power at each time step are defined

RegG(3).Gen_Qref=Gen_Qref(3);
RegG(3).Gen_Pref=Gen_Pref(3);

% Definition of Droop Parameters

RegG(3).V = [0.9,0.97,1.03,1.10]*400;
RegG(3).w = 2*pi*50*[0.98,0.994,1.006,1.02];

RegG(3).Qmax = Gen_Qmax(3);
RegG(3).Qmin = Gen_Qmin(3);

RegG(3).Pmax = Gen_Pmax(3);
RegG(3).Pmin = Gen_Pmin(3);

    % Droop Parameters
RegG(3).mv = -(  (RegG(3).Pmax)    -   (RegG(3).Pmin) )   /  ((RegG(3).V(4)-RegG(3).V(1)));        
RegG(3).mf =  (  (RegG(3).Qmax     -    RegG(3).Qmin) )   /  ((RegG(3).w(4)-RegG(3).w(1)));

%% GENERATOR 4, NODE 11, NATURAL GAS COGENERATOR

RegG(4).nG=11;                                                             % Generator Node

% Definition of the minimum power factor in order to define the
% problem's constraints

RegG(4).cosphimin=0.2;                                                     
RegG(4).b=tan(acos(RegG(4).cosphimin));

% The reference active and reactive power at each time step are defined

RegG(4).Gen_Qref=Gen_Qref(4);
RegG(4).Gen_Pref=Gen_Pref(4);

% Definition of Droop Parameters

RegG(4).V = [0.9,0.97,1.03,1.10]*400;
RegG(4).w = 2*pi*50*[0.98,0.994,1.006,1.02];

RegG(4).Qmax = Gen_Qmax(4);
RegG(4).Qmin = Gen_Qmin(4);

RegG(4).Pmax = Gen_Pmax(4);
RegG(4).Pmin = Gen_Pmin(4);

    % Droop Parameters

RegG(4).mv = -(  ( RegG(4).Pmax    -     RegG(4).Pmin ) )   /  ((RegG(4).V(4)-RegG(4).V(1)));        
RegG(4).mf =  (  ( RegG(4).Qmax     -    RegG(4).Qmin ) )   /  ((RegG(4).w(4)-RegG(4).w(1)));

%% GENERATOR 5, NODE 12, SOLAR PANEL

RegG(5).nG=12;

RegG(5).Gen_Pref=Gen_Pref(5);
RegG(5).Gen_Qref=Gen_Qref(5);
RegG(5).Qmax = Gen_Qmax(5);
RegG(5).Qmin = Gen_Qmin(5);

RegG(5).V =   [1.05,1.15]*400;
RegG(5).w = 2*pi*50*[0.98,0.994,1.006,1.02];

% The droop parameters of the renewable sources are defined inside the
% script regolazione_generatore

%% GENERATOR 6, NODE 13, WIND TURBINE

RegG(6).nG=13;

RegG(6).Gen_Pref=Gen_Pref(6);
RegG(6).Gen_Qref=Gen_Qref(6);

RegG(6).Qmax = Gen_Qmax(6);
RegG(6).Qmin = Gen_Qmin(6);

RegG(6).V =   [1.05,1.15]*400;
RegG(6).w = 2*pi*50*[0.98,0.994,1.006,1.02];

% The droop parameters of the renewable sources are defined inside the
% script regolazione_generatore

%% LOAD 1, NODE 2, RLC LOAD

RegL(1).nL     =      2;
RegL(1).Pr     =    Pr1;                                                   % Active Absorbed Power
RegL(1).Qind   =  Qind1;                                                   % Inductive Reactive Absorbed Power
RegL(1).Qcap   =  Qcap1;                                                   % Capacitive Reactive Absorbed Power

%% LOAD 2, NODE 5, COSTANT LOAD

RegL(2).nL=5;
RegL(2).Sn=SL2;

%% LOAD 3, NODE 6, RLC LOAD

RegL(3).nL   =     6;
RegL(3).Pr   =   Pr2;
RegL(3).Qind = Qind2;
RegL(3).Qcap = Qcap2;

%% LOAD 4, NODE 9, COSTANT LOAD
RegL(4).nL=9;
RegL(4).Sn=SL4;

%% LOAD 5, NODE 11, COSTANT LOAD

RegL(5).nL=11;
RegL(5).Sn=SL5;
