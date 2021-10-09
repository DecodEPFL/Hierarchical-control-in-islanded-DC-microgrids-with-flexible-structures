
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
t_init = 540;
Tsim   = 60;


for t = t_init : (t_init+Tsim-1)
    t
   
 
%         if or(mod(t,MPCpar.t_s) == 0, t == t_start)
%             
%             PL_MPC             =  PL_nom_s(:,  ceil(t/MPCpar.t_s) : ceil(t/MPCpar.t_s )+ MPCpar.N -1);
%   
% %%% Uncomment if you want the MPC to know the measure the real disturbance
% 
%               PL_MPC(:,1)        =  P_nom(:,t);
%           
%             disp('*** MPC started  ***' );
%               [P_ref_MPC, delta_g_k, delta_b_k, Pred_PV_k, Delta_P_slack_max, Delta_P_slack_min,sol_k]  =  HMPC_control(SOC, SOC_0, delta_g_k, delta_b_k, PL_MPC, MPCpar);              
%             disp('*** MPC finished ***' );   
%         end

SOC       =  datalog.SOC(:,t);   
P_ref_MPC =  datalog.MPC.P_ref(1:end-1,t);
delta_g_k =  [P_ref_MPC(1)>0;    P_ref_MPC(4)>0];
Pred_PV_k =  datalog.MPC.Pred_PV(:,t);
Delta_P_slack_max =  datalog.OPF.DeltaPmax(1:end-1,t);
Delta_P_slack_min =  datalog.OPF.DeltaPmin(1:end-1,t);

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



%     if or(mod(t,OPFpar.t_s)==0 , t == t_start)
%             disp('***OPF started***');
%                   OPF;
%             disp('***OPF finished***');
%     end

%Vref = datalog.OPF.Vref(:,t);
v_opf = datalog.OPF.v_opf(:,t);

%     
if  (Pred_PV_k) < 0        
    Vref = [v_opf([sourcenodes_mod pvnodes])];
else
    Vref = [v_opf(sourcenodes_mod); 0];
end 

    if any(abs(Delta_P_slack_opf)>2)
        keyboard
    end
    
% datalog.MPC.sol(:,t-t_start+1)     =   sol_k;

% datalog.MPC.P_ref(:,t-t_start+1)   =   [P_ref_MPC; -P_L(1,1)];

% datalog.MPC.Pred_PV(:,t-t_start+1) =   Pred_PV_k;

% datalog.OPF.v_opf(:,t-t_start+1)   =   v_opf;
% datalog.OPF.Vref (:,t-t_start+1)   =   Vref;
%  
% datalog.OPF.DeltaP(:,t-t_start+1)    = Delta_P_slack_opf;
% datalog.OPF.DeltaPmax(:,t-t_start+1) = Delta_P_slack_max_mod;
% datalog.OPF.DeltaPmin(:,t-t_start+1) = Delta_P_slack_min_mod;
% 
% datalog.OPF.Pg(:,t-t_start+1)        = Delta_P_slack_opf + [P_ref_MPC; -P_L(1,1)];
 

% OPF_sim;
% 
% datalog.OPFsim.v_sim(:,t-t_start+1)      = v_sim_opt;   
% datalog.OPFsim.p_G_sim(:,t-t_start+1)    = p_G_sim_opt;   


if or(and(t>553,t<557),and(t>567,t<573))
    disp('***Simulation started***');
        tic
           networksimulation
        toc
    disp('***Simulation finished***');

    y_plot(:,(t-t_start)*(length(tr))+1 :(t-t_start+1)*(length(tr)) )     =   yout(1:end,:)';
    t_plot(:,(t-t_start)*(length(tr))+1 :(t-t_start+1)*(length(tr)) )     =   time(1:end,:)' + (t-1)*dt_sim ;
end 
%     
%%%%%   P_G_tot_ref_plot =   [P_G_tot_ref_plot, (Delta_P_slack_JSP + p_G_ref)*ones(1,length(time(1:end-1,:)))];

%P_batt =  datalog.OPF.Pg(batt_nodes,t-t_start+1);
%P_batt =   p_G_sim_opt(batt_nodes);

%SOC = SOC - (dt_sim/3600)./(MPCpar.Cb').*( (diag(1./n_dh)*max( P_batt, zeros(nb,1))) +  diag(n_dh)*min( P_batt, zeros(nb,1)) );

    
end

%% Plotting results 

sim_init = 554;
sim_centr = 555;
sim_end  = 556;

sourcenodes_mod = 2:5;
figure
plot(  t_plot(:,(sim_init)*(length(tr))+1  :(sim_end)*(length(tr)))/3600 ,   y_plot(sourcenodes_mod,(sim_init-1)*(length(tr))+1 :(sim_end-1)*(length(tr)) ),'linewidth',2.5)
hold on
plot(  t_plot(:,(sim_init)*(length(tr))+1  :(sim_end)*(length(tr)))/3600 ,  kron(datalog.OPF.v_opf(sourcenodes_mod,sim_init:sim_end-1),ones(1,length(tr))),':k','linewidth',1.5)
xlim( [sim_centr/60-2/(60*60), sim_centr/60+2/(60*60)])
% ylim([min(datalog.OPF.v_opf(sourcenodes,sim_init:sim_end-1),[],1)-0.5, max(datalog.OPF.v_opf(sourcenodes,sim_init:sim_end-1),[],1)+0.5])
grid on


figure
plot(  t_plot(:,(sim_init)*(length(tr))+1  :(sim_end)*(length(tr)))/3600 ,   y_plot(loadnodes,(sim_init-1)*(length(tr))+1 :(sim_end-1)*(length(tr)) ),'linewidth',2.5)
hold on
%plot(  t_plot(:,(569)*(length(tr))+1  :(571)*(length(tr)))/3600 ,  kron(datalog.OPF.v_opf(2,569:570),ones(1,length(tr))),'--k','linewidth',2.5)
xlim( [sim_centr/60-2/(60*60), sim_centr/60+2/(60*60)])
%ylim([min(datalog.OPF.v_opf(2,sim_init:sim_end-1))-0.5, max(datalog.OPF.v_opf(2,sim_init:sim_end-1))+0.5])
grid on
