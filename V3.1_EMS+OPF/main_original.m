
Initialization


%% Simulation 

SOC_0                   =   [0.6;0.5;0.4];
SOC                     =    SOC_0;
delta_g_k               =    zeros(MPCpar.ng,1);
delta_b_k               =    zeros(MPCpar.nb,1);
Pred_PV_k               =    0;
delta_nodes             =    ones(nn,1);
delta_lines             =    ones(mm,1);
n_ch                    =    MPCpar.n_ch;
n_dh                    =    MPCpar.n_dh;
t_init = t_start;
v_opf = v_nom;
% 
% t_init = 9*60;
% Tsim   = 60;
% SOC = datalog.SOC(:,t_init);


for t = t_init : (t_init+Tsim-1)
    t
   
 
datalog.SOC(:,t-t_start+1) = SOC;

        if or(mod(t,MPCpar.t_s) == 0, t == t_start)
            
            PL_MPC             =  PL_nom_s(:,  ceil(t/MPCpar.t_s) : ceil(t/MPCpar.t_s )+ MPCpar.N -1);
  
%%% Uncomment if you want the MPC to know the measure the real disturbance

          %    PL_MPC(:,1)        =  P_nom(:,t);
            tic
            disp('*** MPC started  ***' );
              [P_ref_MPC, delta_g_k, delta_b_k, Pred_PV_k, Delta_P_slack_max, Delta_P_slack_min,sol_k]  =  HMPC_control(SOC, SOC_0, delta_g_k, delta_b_k, PL_MPC, MPCpar);              
            disp('*** MPC finished ***' );   
            t_MPC = toc
        end
    
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

Y_L      =  diag(yl);
I_L      =  diag(il);
P_L      =  diag(pl);
P_L(1,1) =  P_L(1,1) - Pred_PV_k; 

% If the generator is disconnected, the boolean variable of the
% corresponding line is set to 0, otherwise it is set to 1
% 
    if or(mod(t,OPFpar.t_s)==0 , t == t_start)
            disp('***OPF started***');
            tic
                  OPF;
            t_OPF = toc
            disp('***OPF finished***');
    end
%     
if  (Pred_PV_k) < 0        
    Vref = [v_opf([sourcenodes pvnodes])];
else
    Vref = [v_opf(sourcenodes); 0];
end 

%     if any(abs(Delta_P_slack_opf)>2)
%         keyboard
%     end
    
datalog.MPC.sol(:,t-t_start+1)     =   sol_k;

datalog.MPC.P_ref(:,t-t_start+1)   =   [P_ref_MPC; -P_L(1,1)];

datalog.MPC.Pred_PV(:,t-t_start+1) =   Pred_PV_k;

datalog.OPF.v_opf(:,t-t_start+1)   =   v_opf;
datalog.OPF.Vref (:,t-t_start+1)   =   Vref;
 
datalog.OPF.DeltaP(:,t-t_start+1)    = Delta_P_slack_opf;
datalog.OPF.DeltaPmax(:,t-t_start+1) = Delta_P_slack_max_mod;
datalog.OPF.DeltaPmin(:,t-t_start+1) = Delta_P_slack_min_mod;

datalog.OPF.Pg(:,t-t_start+1)        = Delta_P_slack_opf + [P_ref_MPC; -P_L(1,1)];
 


OPF_sim;

datalog.OPFsim.v_sim(:,t-t_start+1)      = v_sim_opt;   
datalog.OPFsim.p_G_sim(:,t-t_start+1)    = p_G_sim_opt;   

% 
% if and(t>566,t<574)
%     disp('***Simulation started***');
%         tic
%             networksimulation
%         toc
%     disp('***Simulation finished***');
% 
%     y_plot(:,(t-t_start)*(length(tr))+1 :(t-t_start+1)*(length(tr)) )     =   yout(1:end,:)';
%     t_plot(:,(t-t_start)*(length(tr))+1 :(t-t_start+1)*(length(tr)) )     =   time(1:end,:)' + (t-1)*dt_sim ;
% end 
%     
%%%%%   P_G_tot_ref_plot =   [P_G_tot_ref_plot, (Delta_P_slack_JSP + p_G_ref)*ones(1,length(time(1:end-1,:)))];

%P_batt =  P_ref_MPC(batt_nodes);
%P_batt =  datalog.OPF.Pg(batt_nodes,t-t_start+1);
P_batt =   p_G_sim_opt(batt_nodes);

SOC = SOC - (dt_sim/3600)./(MPCpar.Cb').*( (diag(1./n_dh)*max( P_batt, zeros(nb,1))) +  diag(n_dh)*min( P_batt, zeros(nb,1)) );

    
end

%% Plotting results 

%load('data_06-12-2020')
plotmain;
