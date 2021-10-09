function     [P_ref, delta_g_k, delta_b_k, Pred_PV_k, Delta_P_slack_max, Delta_P_slack_min, sol_k]  =  MPC_control(SOC_k, SOC_0, delta_g_k, delta_b_k, PL_MPC, MPCpar);

yalmip('clear');

nb = MPCpar.nb;
ng = MPCpar.ng;

nl = size(PL_MPC,1);

N         =   MPCpar.N;
%%
batt_nodes = MPCpar.batt_nodes;
gen_nodes  = MPCpar.gen_nodes;

c_b = MPCpar.c_b;
c_g = MPCpar.c_g;

Cb        =   MPCpar.Cb;
Socmax    =   MPCpar.Socmax;
Socmin    =   MPCpar.Socmin;
SOC_nom   =   MPCpar.SOC_nom;
Pmax      =   MPCpar.Pmax;
Pmin      =   MPCpar.Pmin;
n_ch      =   MPCpar.n_ch;
n_dh      =   MPCpar.n_dh;
t_s       =   MPCpar.t_s;


SOC_N     =   sdpvar(nb,N+1);

P_b_ch    =   sdpvar(nb,N);
P_b_dh    =   sdpvar(nb,N);
P_b       =   sdpvar(nb,N);

delta_b   =   binvar(nb*N,1);

P_gen     =   sdpvar(ng,N);
delta_g   =   binvar(ng*N,1);

%delta_pv   =   binvar(1,N);

Pred_PV  =   sdpvar(1,N);

eps_soc  = sdpvar(nb,1);

Delta_SOC = sdpvar(nb,N);

A         =  ( eye(nb) ) ;

B         =  diag(-t_s/60./Cb');


[Ac, Buc] =   MPC_matmaker(A,B,N); 

%%
J      =          (c_g.*P_gen(:))'     *    (c_g.*P_gen(:))          + ...
       +          (c_b.*P_b(:))'        *   (c_b.*P_b(:))            + ...
       +        5*1e4  *  (Pred_PV)   *   (Pred_PV)'          + ...
       +          1e6  *  (delta_g  -  [delta_g_k; delta_g(1:end-ng)] )'  *  (delta_g -  [delta_g_k; delta_g(1:end-ng)] ) + ...
       +        5*1e7  *  (delta_b  -  [delta_b_k; delta_b(1:end-nb)] )'  *  (delta_b -  [delta_b_k; delta_b(1:end-nb)] ) + ... 
       +      2.5*1e9  *  (eps_soc')       *    eps_soc; % + ...
     %  +           1e3  *  (Delta_SOC(:)') * Delta_SOC(:) ;

vinc = [    SOC_N(:)            ==     Ac*SOC_k  +  Buc*(  P_b(:) ) ;     ...
            Delta_SOC           ==     SOC_N(:,2:end) - kron(ones(1,N), SOC_0);    ...
            P_b                 ==     diag(1./n_dh)*P_b_dh -  diag(n_ch)*P_b_ch; ...
            SOC_N(:,2:end)      <=     kron(Socmax',ones((1),N)) ; ...
            SOC_N(:,2:end)      >=     kron(Socmin',ones((1),N)) ; ... 
            P_b_dh(:)           <=     delta_b.*kron(ones(N,1),Pmax(batt_nodes)'); ... 
            P_b_ch(:)           <=    -(ones(nb*N,1) - delta_b).*kron(ones(N,1),Pmin(batt_nodes)'); ...
            P_b_dh(:)           >=     zeros(nb*N,1); ... 
            P_b_ch(:)           >=     zeros(nb*N,1); ... 
            P_gen(:)            <=     delta_g .* kron(ones(N,1), Pmax(gen_nodes)'); ... 
            P_gen(:)            >=     delta_g .* kron(ones(N,1), Pmin(gen_nodes)'); 
sum(reshape(delta_g,[ng,N]),1)  >=     ones(1,N); ...
            Pred_PV             <=     zeros(1,N); ...
            Pred_PV             >=     PL_MPC(1,1:N); ...
            SOC_N(:,end)        <=     SOC_0  + eps_soc; ...
            SOC_N(:,end)        >=     SOC_0  - eps_soc; ...
            eps_soc             >=     zeros(nb,1); ...
            zeros(1,N)          ==     ones(1,nb)*(P_b_dh - P_b_ch) + ones(1,ng)*P_gen - ones(1,nl)*PL_MPC(:,1:N) + Pred_PV; ...
       ];
               
opt = sdpsettings('solver','cplex','showprogress',0,'verbose',0);
opt.cplex.timelimit = 360;
tstart = tic;
sol = solvesdp(vinc,J,opt)
tstop = toc(tstart);
%
%  if tstop > opt.cplex.timelimit
%     keyboard
% opt = sdpsettings('solver','bnb','showprogress',0,'verbose',0);
% tstart2 = tic;
% sol = solvesdp(vinc,J,opt)
% tstop2 = toc(tstart2);
%     if tstop2 >120
%      keyboard
%     end
%  end

if sol.problem >0
    if (sol.problem ~=3)
    keyboard;
    end
end
sol_k = sol.problem;

P_ref    = zeros(nb+ng,1);

P_ref(MPCpar.batt_nodes)    = [ double(P_b_dh(:,1)) - double(P_b_ch(:,1))];
P_ref(MPCpar.gen_nodes)     = [ double(P_gen(:,1))];

delta_g_k = double(delta_g(1:ng));
delta_b_k = double(delta_b(1:nb));

Pred_PV_k = double(Pred_PV(1));

% Definition of the bounds for the slack variables that the OPF can
% manipulate

P_b_max   = min ( Pmax(MPCpar.batt_nodes)', (n_dh).*Cb'.*(double(SOC_N(:,1)) -  Socmin')./(t_s/60));
P_b_min   = max ( Pmin(MPCpar.batt_nodes)', (1./n_ch).*Cb'.*(double(SOC_N(:,1))  -  Socmax')./(t_s/60));
P_g_max   = Pmax( MPCpar.gen_nodes)';
P_g_min   = Pmin( MPCpar.gen_nodes  )';

Delta_P_slack_max = zeros(nb+ng,1);
Delta_P_slack_min = zeros(nb+ng,1);

Delta_P_slack_max(MPCpar.batt_nodes) = P_b_max - P_ref(MPCpar.batt_nodes);
Delta_P_slack_min(MPCpar.batt_nodes) = P_b_min - P_ref(MPCpar.batt_nodes);

Delta_P_slack_max(MPCpar.gen_nodes)  = P_g_max.*delta_g_k - P_ref(MPCpar.gen_nodes);
Delta_P_slack_min(MPCpar.gen_nodes)  = P_g_min.*delta_g_k - P_ref(MPCpar.gen_nodes);

end