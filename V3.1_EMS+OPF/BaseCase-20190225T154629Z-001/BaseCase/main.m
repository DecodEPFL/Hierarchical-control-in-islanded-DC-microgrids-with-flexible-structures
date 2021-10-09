%% Initialization of Simulation

Initialization;

prompt1=('Do you want to start a NEW simulation (1) or to continue an OLD one (2) ? ');
a1=input(prompt1);                                                 

if (a1~=1 && a1~=2)
    
display('Error!!!');

elseif(a1==2)

    load datalog
 
    prompt2=('Initial time step to simulate ( 1 time step = 30 seconds) ?   ');
    a=input(prompt2);
    
    prompt3=('Final time step to simulate   ?   ');
    T=input(prompt3);
    
    % The data concerning the last step of the "OLD simulation"
    
        Soc = datalog.X(end-nb+1 : end, a);
        v = datalog.X(end-nb, a);

        Gen_Pref = [datalog.U(1:4,a)', real(Pgen(a,5)), real(Pgen(a,6))];
        Gen_Qref = [datalog.U(5:8,a)', imag(Pgen(a,5)), imag(Pgen(a,6))];
        
        Gen_Pmax = datalog.Gen_Pmax(:,a)';
        Gen_Pmin = datalog.Gen_Pmin(:,a)';

        Pr1=Pr_profile(a,1);
        Pr2=Pr_profile(a,2);

        Qind1=Qind_profile(a,1);
        Qind2=Qind_profile(a,2);

        Qcap1=Qcap_profile(a,1);
        Qcap2=Qcap_profile(a,2);

        SL2=Pload(a,2);
        SL4=Pload(a,4);
        SL5=Pload(a,5);

        generatori_e_carichi;
        
       
elseif (a1==1)
    
    prompt2=('Initial time step to simulate ( 1 time step = 30 seconds) ?   ');
    a=input(prompt2);
    
    prompt3=('Final time step to simulate  ?   ');
    T=input(prompt3);   
end

% Some variables are however initialized to zero
Umpc = zeros(8,N);
hsoft=0;
Jopt=0;
Socmpc=Soc;
PGEN = zeros(8,1);

%% Simulation

% I simulate the controlled system at each istant apart from the last time 
% one where the MPC is not performed

for t = a : (T-1)         
    
    t
    
    % I set simulator options
    options.IterMax=100;
    options.IterOk=2;
    options.Tol=[1e-5,1e-5,0.10];
    
    %%% the network is simulated with NR
    [results,details] = island_PF(netdata,RegG,RegL,options);
    close figure 100
    
    dv=details.dv(end);
    dw=details.dw(end);
    ds=details.ds(end);
  
    % convergence of simulator and number of iterations
    conv = details.convergence;
    iter = details.iterations;
   
    % If the simulator does not converge the simulation is interrupted
    
    if (conv==0)
        load handel
        player = audioplayer(y, Fs);
        play(player);
        keyboard
    end
    
  % The voltages, phases and frequency are saved
  
    V           =  abs(results.V);
    angleV      =  angle(results.V);
    w           =  results.w;
    
    
    %%% Through the active generated power at that istant the actual Soc of each battery is computed. 
      % Morover, the maximum and minimum powers that batteries can deliver are also computed.
    
    [Gen_Pmax, Gen_Pmin, Soc] = batt_model( Soc, real(details.Sgen), gen_batt, k, MPCpar, Gen_Pmax, Gen_Pmin, Gen_Pmax0, Gen_Pmin0, nodi_ctrl_gen);
        
    
    %%% Definition of the state vector
    x = [V; angleV(2:end); w; v];   
   
    
    %%% Definition of the actual active and reactive generated powers of
    %%% controllable sources
    
    P0 = [real(details.Sgen(1:4));imag(details.Sgen(1:4))];

    %%% Now, the network model is linearized with the actual state and so the
    %%% jacobian can be computed

    JinvReal = inv(RealJacobian(V,angleV,w,RegG,RegL,netdata));            % Jacobian based on the real network paramenters
    % rcond(RealJacobian(V,angleV,w,RegG,RegL,netdata))


  % The MPC is set to be performed every "s" steps
    s = 2;                                                                 % s=2 means to run the MPC at each minute
   

if (mod(t,s)==0)
   
            %%% MPC CONTROL
             [Umpc, v, Jopt, val_h, TotU, Soc, PGEN, hsoft, exitflag] = MPC_control_qp(x, Soc, JinvReal, var_D(:,(t+1):(t+N)), netdata, MPCpar, Gen_Pref(1:4), Gen_Qref(1:4), k, gen_batt, Gen_Pmax0, Gen_Pmin0, Gen_Qmax0(1:4), Gen_Qmin0(1:4), RegG, P0);
   
             
            %%% Definition of reference power after the MPC control action
            
            Gen_Pref(1)     =   Umpc(1,1)                                      +        Gen_Pref(1);

            Gen_Qref(1)     =   Umpc(1     + length(nodi_ctrl_gen),1)          +        Gen_Qref(1);

            Gen_Pref(2)     =   Umpc(2,1)                                      +        Gen_Pref(2);

            Gen_Qref(2)     =   Umpc(2     + length(nodi_ctrl_gen),1)          +        Gen_Qref(2);

            Gen_Pref(3)     =   Umpc(3,1)                                      +        Gen_Pref(3);

            Gen_Qref(3)     =   Umpc(3     + length(nodi_ctrl_gen),1)          +        Gen_Qref(3);

            Gen_Pref(4)     =   Umpc(4,1)                                      +        Gen_Pref(4);

            Gen_Qref(4)     =   Umpc(4     + length(nodi_ctrl_gen),1)          +        Gen_Qref(4);    
    
end

    %%% Definition reference powers of uncontrollable sources for the next
    %%% time step to simulate
    
    Gen_Pref(5)=real(Pgen(t+1,5));
    Gen_Qref(5)=imag(Pgen(t+1,5));
 
    Gen_Pref(6)=real(Pgen(t+1,6));
    Gen_Qref(6)=imag(Pgen(t+1,6));
    
 	%%% Definition of loads powers for the next time step to simulate
    
    Pr1=Pr_profile(t+1,1);
    Pr2=Pr_profile(t+1,2);
    
    Qind1=Qind_profile(t+1,1);
    Qind2=Qind_profile(t+1,2);
    
    Qcap1=Qcap_profile(t+1,1);
    Qcap2=Qcap_profile(t+1,2);
 
    SL2=Pload(t+1,2);
    SL4=Pload(t+1,4);
    SL5=Pload(t+1,5);
   
    generatori_e_carichi;

    %%% Data saving
    
    datalog.U(1:4,t+1)        =   (Gen_Pref(:,1:4))';
    datalog.U(5:8,t+1)        =   (Gen_Qref(:,1:4))';
    
    datalog.DU(:,t)           =   Umpc(:,1);
    datalog.hsoft(:,t)        =   hsoft;
    datalog.Jopt(1,t)         =   Jopt;

    datalog.Pgen(:,t+1)       =   PGEN;
    datalog.Socmpc(:,t+1)     =   Socmpc;
          
    datalog.Gen_Pmax(:,t+1)   =   Gen_Pmax';
    datalog.Gen_Pmin(:,t+1)   =   Gen_Pmin';
   
    datalog.X(:,t)            =   [x; Soc];
    datalog.Sgen(:,t)         =   (details.Sgen);
    datalog.Sload(:,t)        =   (details.Sload);
    

    save ('datalog','datalog')

end

% Finally the last time step is simulated without the MPC control action

if(t==(T-1))
    
        t=t+1          
       
        
   [results,details]=island_PF(netdata,RegG,RegL,options);
   close figure 100
   
            V           =  abs(results.V);
            angleV      =  angle(results.V);
            w           =  results.w;

            x=[V;angleV(2:end);w;v];                                
            
    [Gen_Pmax, Gen_Pmin, Soc] = batt_model( Soc, real(details.Sgen), gen_batt, k, MPCpar, Gen_Pmax, Gen_Pmin,Gen_Pmax0, Gen_Pmin0, nodi_ctrl_gen);
    
    datalog.X(:,t)     = [x;Soc];
    datalog.Sgen(:,t)  = details.Sgen;
    datalog.Sload(:,t) = details.Sload;
            
          
end

%% End of the Simulation       
       
save ('datalog','datalog')
