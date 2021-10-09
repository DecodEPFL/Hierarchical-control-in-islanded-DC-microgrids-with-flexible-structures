function [Umpc, Jopt, h, SOC]=MPC_control(x, Jinv, Jloss, Ploss, var_D, netdata, MPCpar, Umpc_0,t)

 
n_nodes          =   netdata.n_nodes;
nodi_ctrl_gen    =   netdata.nodi_ctrl_gen;
nodi_batt_pq     =   netdata.nodi_batt_pq;
n_b              =   length(nodi_batt_pq);

nu = length(nodi_ctrl_gen)*2;                                                % Numero ingressi
nx = length(x);                                                              % Numero stati
nd = (n_nodes-1)*2;            %Tolgo lo slack
% 
% ro          = MPCpar.Cf.ro;
% p_int       = MPCpar.Cf.p_int;
% p_delta     = MPCpar.Cf.p_delta;

% p_u         = MPCpar.Cf.p_u;
% 
N           = MPCpar.N;
% 
Vmax        = MPCpar.Vmax;
Vmin        = MPCpar.Vmin;
% Umax        = MPCpar.Umax;
% Umin        = MPCpar.Umin;

%% Creazione matrici
%
A = eye(nx);
% Aggiungo alla matrice A la dinamica del SOC delle batterie
A(2*n_nodes + 2*n_b + 1, 2*n_nodes + 1)  =  - MPCpar.step/MPCpar.C_gen1;
A(2*n_nodes + 2*n_b + 2, 2*n_nodes + 2)  =  - MPCpar.step/MPCpar.C_gen2;
% Aggiungo la dinamica del SOC slack
A(end, 2*n_nodes-1)  =  - MPCpar.step/MPCpar.C_ns;
%
%
B=Jinv(:, [(nodi_ctrl_gen-1), ((nodi_ctrl_gen-1) + (n_nodes-1))] );
% Aggiungo la dinamica delle potenze di slack che dipendono dal bilancio di
% potenza, quindi sommano la variazione di potenza dei generatori, sia attiva che reattiva
B = [B; -1*ones(1,length(nodi_ctrl_gen)), zeros(1,length(nodi_ctrl_gen)); zeros(1,length(nodi_ctrl_gen)), -1*ones(1,length(nodi_ctrl_gen))];
%
% Aggiungo le righe per le potenze batterie attive
B = [B; zeros(n_b,size(B,2))];
B(end-1,1) = 1;
B(end,2)   = 1;
% Aggiungo le righe per le potenze batterie reattive
B = [B; zeros(n_b,size(B,2))];
B(end-1,3) = 1;
B(end,4)   = 1;
% Aggiungo le righe per i SOC_pq
B = [B; zeros(n_b,size(B,2))];
B(end-1,1) = - MPCpar.step/MPCpar.C_gen1;
B(end,2)   = - MPCpar.step/MPCpar.C_gen2;
% Aggiungo dinamica del SOC_s
B = [B; +MPCpar.step/MPCpar.C_ns*ones(1,length(nodi_ctrl_gen)), zeros(1,length(nodi_ctrl_gen))];
%
M=Jinv;
% Aggiungo dinimica relativa alla potenza di slack
M = [M; ones(1,n_nodes-1), zeros(1,n_nodes-1); zeros(1,n_nodes-1), ones(1,n_nodes-1)];
% Aggiungo righe di zeri per quanti riguarda le potenze delle
% batterie PQ
M = [M; zeros(2*n_b,size(M,2))];
%Aggiungo dinamica SOC_pq
M = [M; zeros(n_b,size(M,2))];
%Aggiungo dinamica SOC_s
M = [M; -MPCpar.step/MPCpar.C_ns*ones(1,n_nodes-1), zeros(1,n_nodes-1)];
%
%
%
% Aggiungo modellistica perdite
A = [A, zeros(size(A,1),n_nodes-1); zeros(n_nodes-1,size(A,2)), eye(n_nodes-1)];
B = [B; Jloss*Jinv(:, [(nodi_ctrl_gen-1), ((nodi_ctrl_gen-1) + (n_nodes-1))] )];
M = [M; Jloss*Jinv];

[Ac, Buc, Mdc] = MPC_matmaker(A,B,M,N); 
%% Optimization problem


x0          =    sdpvar(nx+length(Ploss),1,'full');
d_U         =    sdpvar(nu,N,'full');
d_D         =    sdpvar(nd,N,'full');
h_v         =    sdpvar(n_nodes-1,N,'full');
h_v_stp     =    sdpvar(n_nodes-1,N,'full');
h_s         =    sdpvar(1,N,'full');
h_u_up      =    sdpvar(1,n_b*2*N,'full');
h_u_dwn     =    sdpvar(1,n_b*2*N,'full');
h_soc_slack =    sdpvar(1,N,'full'); 
X            =    reshape(Ac*x0 + Buc*d_U(:) + Mdc*d_D(:), nx+length(Ploss),  []   );

X_V               =    X(1:n_nodes-1,2:end);
X_Ps              =    X(2*n_nodes-1,2:end);
X_Qs              =    X(2*n_nodes,2:end);
X_Pbatt_1         =    X(2*n_nodes+1,2:end);
X_Pbatt_2         =    X(2*n_nodes+2,2:end);
X_Qbatt_1         =    X(2*n_nodes+n_b+1,2:end);
X_Qbatt_2         =    X(2*n_nodes+n_b+2,2:end);
X_SOC_pq_1        =    X(2*n_nodes+2*n_b+1,2:end);
X_SOC_pq_2        =    X(2*n_nodes+2*n_b+2,2:end);
X_SOC_s           =    X(2*n_nodes+2*n_b+3,2:end);
X_Ploss           =    X(2*n_nodes+2*n_b+4:end,1:end);
% 
% var_X_V_2         =    X_V(2,:) - netdata.Vnom*ones(1,N);
% var_X_V_3         =    X_V(3,:) - netdata.Vnom*ones(1,N);
% var_X_V_4         =    X_V(4,:) - netdata.Vnom*ones(1,N);
% var_X_V_5         =    X_V(5,:) - netdata.Vnom*ones(1,N);



J_V    = 1e2*MPCpar.Cf.p_v(2)*((h_v_stp(:)')*(h_v_stp(:))); % ( MPCpar.Cf.p_v(2).*(var_X_V_2*var_X_V_2') + MPCpar.Cf.p_v(3).*(var_X_V_3*var_X_V_3') )+ MPCpar.Cf.p_v(4)*(var_X_V_4*var_X_V_4') +  MPCpar.Cf.p_v(5)*(var_X_V_5*var_X_V_5');
J_h    =     MPCpar.Cf.h*(sum(h_v(:)) + 1e2*sum(h_soc_slack(:)) + sum(h_s) + sum(h_u_up) + sum(h_u_dwn));
J_dU   =     MPCpar.Cf.dU(1)*d_U(1,:)*d_U(1,:)'+ MPCpar.Cf.dU(2)*d_U(2,:)*d_U(2,:)' + MPCpar.Cf.dU(3)*d_U(3,:)*d_U(3,:)'+ MPCpar.Cf.dU(4)*d_U(4,:)*d_U(4,:)';
J_loss =     MPCpar.Cf.Ploss*(X_Ploss(:)'*X_Ploss(:));
% J_P    =  (  MPCpar.Cf.P(1)*X_Pbatt_1(:)'*X_Pbatt_1(:) + MPCpar.Cf.P(2)*X_Pbatt_2(:)'*X_Pbatt_2(:) + MPCpar.Cf.Q(1)*X_Qbatt_1(:)'*X_Qbatt_1(:) + MPCpar.Cf.Q(2)*X_Qbatt_2(:)'*X_Qbatt_2(:));% + MPCpar.Cf.P(3)*X_Ps(:)'*X_Ps(:)+ MPCpar.Cf.Q(3)*X_Qs(:)'*X_Qs(:));

% J_soc  =    1e1*(X_SOC_s(2:end) - 0.5*ones(1,N))*(X_SOC_s(2:end) - 0.5*ones(1,N))';% + ...
%       1e1*(X_SOC_pq_1(1:end) - 0.5*ones(1,N))*(X_SOC_pq_1(1:end) - 0.5*ones(1,N))' + ...
%       1e1*(X_SOC_pq_2(1:end) - 0.5*ones(1,N))*(X_SOC_pq_2(1:end) - 0.5*ones(1,N))' ;

% if t > 1
% J_var  = 1e2*MPCpar.Cf.dU(1).*(d_U(1,1) - Umpc_0(1))^2 +  1e2*MPCpar.Cf.dU(2).*(d_U(2,1) - Umpc_0(2))^2 +  1e1*MPCpar.Cf.dU(3).*(d_U(3,1) - Umpc_0(3))^2 +  1e1*MPCpar.Cf.dU(3).*(d_U(4,1) - Umpc_0(4))^2;
% else 
% J_var  = 0;
% end

J =     J_V  + ...
        J_h  + ...
        J_dU + ...
      J_loss ;

vinc= [  x0                   ==      [x;Ploss]; ...
         d_D(:)               ==      reshape(var_D,size(d_D(:),1),1);...
         d_U(:)               <=      0.2*ones(n_b*2*N,1) + h_u_up';...
         d_U(:)               >=     -0.2*ones(n_b*2*N,1) - h_u_dwn';...
         X_V(:)               <=      repmat(Vmax,N*(n_nodes-1),1) + h_v(:);...
         X_V(:)               >=      repmat(Vmin,N*(n_nodes-1),1) - h_v(:);...
         X_V(:)               <=      repmat(1.01,N*(n_nodes-1),1) + h_v_stp(:);...
         X_V(:)               >=      repmat(0.99,N*(n_nodes-1),1) - h_v_stp(:);...
         X_Pbatt_1(:)         <=      MPCpar.Pmax_gen1_pu*ones(N,1);...
         X_Pbatt_1(:)         >=      MPCpar.Pmin_gen1_pu*ones(N,1);...
         X_Pbatt_2(:)         <=      MPCpar.Pmax_gen2_pu*ones(N,1);...
         X_Pbatt_2(:)         >=      MPCpar.Pmin_gen2_pu*ones(N,1);...
         X_Qbatt_1(:)         <=      MPCpar.Qmax_gen1_pu*ones(N,1);...
         X_Qbatt_1(:)         >=      MPCpar.Qmin_gen1_pu*ones(N,1);...
         X_Qbatt_2(:)         <=      MPCpar.Qmax_gen2_pu*ones(N,1);...
         X_Qbatt_2(:)         >=      MPCpar.Qmin_gen2_pu*ones(N,1);...
         X_SOC_pq_1(:)        <=      MPCpar.Socmax_gen1*ones(N,1); ...
         X_SOC_pq_1(:)        >=      MPCpar.Socmin_gen1*ones(N,1);...
         X_SOC_pq_2(:)        <=      MPCpar.Socmax_gen2*ones(N,1);...
         X_SOC_pq_2(:)        >=      MPCpar.Socmin_gen2*ones(N,1); ...
         X_SOC_s(:)           <=      MPCpar.Socmax_ns*ones(N,1);...
         X_SOC_s(:)           >=      MPCpar.Socmin_ns*ones(N,1); ...
         X_SOC_s(:)           <=      0.8*ones(N,1) + h_soc_slack';...
         X_SOC_s(:)           >=      0.2*ones(N,1) - h_soc_slack';...
         X_Ps(:)              <=     (MPCpar.Pmax_ns_pu - MPCpar.MaxPLoadStep_up_pu)*ones(N,1)  + h_s'; ...
         X_Ps(:)              >=     (MPCpar.Pmin_ns_pu + MPCpar.MaxPLoadStep_dwn_pu)*ones(N,1) - h_s'; ...
         X_Qs(:)              <=     (MPCpar.Qmax_ns_pu - MPCpar.MaxQLoadStep_up_pu)*ones(N,1)  + h_s'; ...
         X_Qs(:)              >=     (MPCpar.Qmin_ns_pu + MPCpar.MaxQLoadStep_dwn_pu)*ones(N,1) - h_s'; ...
         X_Ps(:)              <=     (X_SOC_s(1:end)' - MPCpar.Socmin_ns*ones(N,1))*MPCpar.C_ns/MPCpar.step - (MPCpar.MaxPLoadStep_up_pu)*ones(N,1)+ h_s'; ...
         X_Ps(:)              >=     (X_SOC_s(1:end)' - MPCpar.Socmax_ns*ones(N,1))*MPCpar.C_ns/MPCpar.step + (MPCpar.MaxPLoadStep_dwn_pu)*ones(N,1)- h_s'; ...         
         h_v(:)               >=      zeros(size(h_v(:),1),1); ...
         h_soc_slack(:)       >=      zeros(size(h_soc_slack(:),1),1);...
         h_s'                 >=      zeros(size(h_s,2),1);...
       ];
%     
% opt=sdpsettings('solver','cplex','showprogress',0,'verbose',0);
opt=sdpsettings('solver','cplex','showprogress',0,'verbose',0);
sol=solvesdp(vinc,J,opt)

if sol.problem >0
    if sol.problem ~=4
    keyboard;
    end
end
%%
Umpc=double(d_U);
Jopt=double(J);
h=[double(h_v);double(h_soc_slack);double(h_s)];

SOC = [double(X_SOC_pq_1(1)); double(X_SOC_pq_2(1)); double(X_SOC_s(1))];
end

 