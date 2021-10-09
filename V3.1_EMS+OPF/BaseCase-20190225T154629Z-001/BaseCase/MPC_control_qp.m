function [Umpc, v, Jopt, val_h, TotU, Soc, Pgen, val_hpsoft, exitflag, V,Aineq_,bineq_,Aeq_,beq_]=MPC_control_qp(x, Soc, Jinv, var_D, netdata, MPCpar, Gen_Pref, Gen_Qref, k, gen_batt,Gen_Pmax,Gen_Pmin,Gen_Qmax, Gen_Qmin, RegG, P0)

%% MPC Initialization

n_nodes          =   netdata.n_nodes;
nodi_ctrl_gen    =   netdata.nodi_ctrl_gen;
wref             =   netdata.wn;    
Vn               =   netdata.Vn;
N                =   MPCpar.N;

ncg              =   length(nodi_ctrl_gen);
nb               =   length(gen_batt);

% The variation with respect to the previous time instant is set as
% disturbance

Dprof = var_D;


% Cost function weights

ro            =   MPCpar.Cf.ro;
rovar         =   MPCpar.Cf.rovar;
rosoc         =   MPCpar.Cf.rosoc;

p_int         =   MPCpar.Cf.p_int;
p_v           =   MPCpar.Cf.p_v;

p_u1_p        =   MPCpar.Cf.p_u1_p;
p_u1_q        =   MPCpar.Cf.p_u1_q;

p_u2_p        =   MPCpar.Cf.p_u2_p;
p_u2_q        =   MPCpar.Cf.p_u2_q;

p_u3_p        =   MPCpar.Cf.p_u3_p;
p_u3_q        =   MPCpar.Cf.p_u3_q;

p_u4_p        =   MPCpar.Cf.p_u4_p;
p_u4_q        =   MPCpar.Cf.p_u4_q;

p_u1_dp       =   MPCpar.Cf.p_u1_dp;
p_u1_dq       =   MPCpar.Cf.p_u1_dq;

p_u2_dp       =   MPCpar.Cf.p_u2_dp;
p_u2_dq       =   MPCpar.Cf.p_u2_dq;

p_u3_dp       =   MPCpar.Cf.p_u3_dp;
p_u3_dq       =   MPCpar.Cf.p_u3_dq;

p_u4_dp       =   MPCpar.Cf.p_u4_dp;
p_u4_dq       =   MPCpar.Cf.p_u4_dq;


% Constraints bounds

Vmax        = MPCpar.Vmax;
Vmin        = MPCpar.Vmin;

SocMax      = MPCpar.SocMax;
SocUB       = MPCpar.SocUB;

SocLB       = MPCpar.SocLB;
SocMin      = MPCpar.SocMin;

uvarnom     = MPCpar.uvarnom;


%% Definition of matrices for Voltages, Phases, Frequency

A=eye(2*n_nodes);                                                                 
B=Jinv(:,[nodi_ctrl_gen, (nodi_ctrl_gen+n_nodes)]);
M=Jinv;

%% Modification of matrices because of FREQUENCY INTEGRATOR

A_tilde1  =  [A , zeros(2*n_nodes,1); zeros(1,2*n_nodes+1)];    
A_tilde1(end,end-1) = -1; A_tilde1(end,end) = 1;

B_tilde1  =  [B;-B(end,:)];

M_tilde1  =  [M zeros(2*n_nodes,1) ; [-M(end,:) 1]];

Dprof_tilde1    =  [Dprof; ones(1,size(var_D,2))*wref];  

nx_tilde1=2*n_nodes+1;

%% Modification of matrices because of REFERENCE POWERS of generators

% The model dynamic to add is Pref(k+1)=Pref(k)+DPref(k)

A_tilde2 = [A_tilde1, zeros(nx_tilde1, 2*ncg); zeros( 2*ncg, nx_tilde1), eye(2*ncg)];
B_tilde2 = [B_tilde1;eye(2*ncg)];
M_tilde2 = [M_tilde1;zeros(2*ncg,nx_tilde1)];

Dprof_tilde2  =  Dprof_tilde1;
nx_tilde2     =  nx_tilde1+2*ncg;

%% Modification of matrices because of ESTIMATED GENERATED POWERS

[A_tilde3, B_tilde3, M_tilde3, nx_tilde3, Dprof_tilde3] = addPowerstates(A_tilde2, B_tilde2, M_tilde2, RegG, N, nodi_ctrl_gen, n_nodes, x, nx_tilde2, nx_tilde1, Dprof_tilde2, Vn);


%% Modification of matrices because of SOCs

[A_tilde4, B_tilde4, M_tilde4, nx_tilde4, Dprof_tilde4] = addSOCstates(A_tilde3, B_tilde3, M_tilde3, gen_batt, nodi_ctrl_gen, nx_tilde3, k , Dprof_tilde3);

%% MPC MATRICES

% I compute the matrices of the system model for the whole prediction
% horizon

[Ac, Buc, Mdc] = MPC_matmaker(A_tilde4, B_tilde4, M_tilde4, N);

% I create a vector to define in which node each battery is placed

pos_batt=zeros(1,nb);

for b=1:nb
    for j=1:ncg
      if (nodi_ctrl_gen(j) == gen_batt(b))
        pos_batt(b) = j;
      end
    end
end


%% Costruction of Aeq e beq

xzero=[x; Gen_Pref'; Gen_Qref'; P0; Soc];
nu=size(Buc,2)/N;

n_h = 3; % Three Slack Variables


%%%%% System Dynamic Constraints
Aeq_ = [eye(size(Ac,1)), -Buc, zeros(size(Ac,1),n_h)];
beq_ = Ac * xzero + Mdc * Dprof_tilde4(:);

%%% Slack Variables' Indeces
i_h=size(Aeq_,2)-n_h+1;
i_hpsoft=size(Aeq_,2)-n_h+2;
i_hsoc=size(Aeq_,2)-n_h+3;


%%%%% VARIABLES' CONSTRAINTS
Aineq_=[];
bineq_=[];

%%%%%%%% VOLTAGES' BOUNDS

for k=1:nx_tilde4:nx_tilde4*(N+1)      
    for t=0:n_nodes-1                
        jj=k+t;
        
        temporaneo=zeros(2,size(Aeq_,2));
        temporaneo(:,jj)=[1 -1]';
        temporaneo(:,i_h)=[-1 -1]';
        
        temporaneo_b=[Vmax(t+1) -Vmin(t+1)]';
        
        Aineq_=[Aineq_;temporaneo];
        bineq_=[bineq_;temporaneo_b];
    end
end

temp_h=zeros(1,size(Aeq_,2));
temp_h(i_h)=-1;
Aineq_=[Aineq_; temp_h];
bineq_=[bineq_;0];

%%%%%%%% SOCs' BOUNDS

for k=(nx_tilde3+1+nx_tilde4):nx_tilde4:nx_tilde4*(N+1)  
  
    for t=0:nx_tilde4-nx_tilde3-1                         
        jj=k+t;
        
        temporaneo=zeros(4,size(Aeq_,2));
        temporaneo(:,jj)=[1 -1 1 -1]';
        temporaneo(:,i_hsoc)=[0 0 -1 -1]';
        
        temporaneo_b=[SocMax(t+1) -SocMin(t+1) SocUB(t+1) -SocLB(t+1)]';
        
        Aineq_=[Aineq_;temporaneo];
        bineq_=[bineq_;temporaneo_b];
    end
end

temp_h=zeros(1,size(Aeq_,2));
temp_h(i_hsoc)=-1;
Aineq_=[Aineq_; temp_h];
bineq_=[bineq_;0];

%%%%%%%% REFERENCE POWERS LIMITS

for k=(nx_tilde1+1+nx_tilde4):nx_tilde4:nx_tilde4*(N+1)
    for t=0:nx_tilde2-nx_tilde1-1
        jj=k+t;
        
        temporaneo=zeros(2,size(Aeq_,2));
        temporaneo(:,jj)=[1 -1]';

        Aineq_=[Aineq_;temporaneo];
    end
end

temporaneo_b=repmat(kron([Gen_Pmax';Gen_Qmax'],[1;0])-kron([Gen_Pmin';Gen_Qmin'],[0;1]),N,1);
bineq_=[bineq_;temporaneo_b];

% An additional constraint is added with regards to the power factor limits
% of the generator number 4

for k=(nx_tilde1+1+nx_tilde4):nx_tilde4:nx_tilde4*(N+1)
   
    t = 7;
    jj  =  k+t;
        
        temporaneo=zeros(2,size(Aeq_,2));
        temporaneo(:,jj)   =  [1 -1]';
        temporaneo(:,jj+4) =  -RegG(4).b*[1 1]'; 

        Aineq_=[Aineq_;temporaneo];
   
end
    
temporaneo_b=repmat([0 0]',N,1);
bineq_=[bineq_;temporaneo_b];
   

%%%%%%%% ESTIMATED GENERATED POWERS LIMITS

for k=(nx_tilde4+nx_tilde2+1):nx_tilde4:nx_tilde4*(N+1)
    for t=0:nx_tilde3-nx_tilde2-1
        jj=k+t;
        
        temporaneo=zeros(2,size(Aeq_,2));
        temporaneo(:,jj)=[1 -1]';

        Aineq_=[Aineq_;temporaneo];
    end
end

temporaneo_b=repmat(kron([Gen_Pmax';Gen_Qmax'],[1;0])-kron([Gen_Pmin';Gen_Qmin'],[0;1]),N,1);
bineq_=[bineq_;temporaneo_b];

% An additional constraint is added with regards to the power factor limits
% of the generator number 4

for k=(nx_tilde2+1+nx_tilde4):nx_tilde4:nx_tilde4*(N+1)
   
    t = 7;
    jj  =  k+t;
        
        temporaneo=zeros(2,size(Aeq_,2));
        temporaneo(:,jj)   =  [1 -1]';
        temporaneo(:,jj-4) =  -RegG(4).b*[1 1]'; 

        Aineq_=[Aineq_;temporaneo];
   
end
    
temporaneo_b=repmat([0 0]',N,1);
bineq_=[bineq_;temporaneo_b];

%%%%%%%% REFERENCE POWER VARIATIONS' BOUNDS

idx_gen=setdiff(1:2*ncg,pos_batt);
for iu=1:numel(idx_gen)
    for k=(nx_tilde4*(N+1)+idx_gen(iu)):nu:(nx_tilde4*(N+1)+nu*N)
        
        temporaneo=zeros(2,size(Aeq_,2));
        temporaneo(:,k)=[1 -1]';

        temporaneo_b=[1 1]'*uvarnom;

        Aineq_=[Aineq_;temporaneo];
        bineq_=[bineq_;temporaneo_b];
    end
end

%%%%%%%% REFERENCE POWER VARIATIONS' BOUNDS FOR BATTERIES

for ib=1:numel(pos_batt)
    for k=(nx_tilde4*(N+1)+pos_batt(ib)):nu:(nx_tilde4*(N+1)+nu*N)
        
        temporaneo=zeros(2,size(Aeq_,2));
        temporaneo(:,k)=[1 -1]';
        temporaneo(:,i_hpsoft)=[-1 -1]';

        temporaneo_b=[1 1]'*uvarnom;

        Aineq_=[Aineq_;temporaneo];
        bineq_=[bineq_;temporaneo_b];
    end
end

temp_h=zeros(1,size(Aeq_,2));
temp_h(i_hpsoft)=-1;
Aineq_=[Aineq_; temp_h];
bineq_=[bineq_;0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COST FUNCTION

f_=zeros(1,size(Aeq_,2));

% Slack Variables' weights
f_(i_h)=ro;
f_(i_hpsoft)=rovar;
f_(i_hsoc)=rosoc;

% Weights
pesiH=zeros(1,size(Aeq_,2));                                

pesiH((2*n_nodes+1):nx_tilde4:nx_tilde4*(N+1))=p_int;                      % Frequency Integrator Error

for iv=1:n_nodes
    pesiH(iv:nx_tilde4:nx_tilde4*(N+1))=p_v;                               % deltaV (errore al quadrato)
    f_(iv:nx_tilde4:nx_tilde4*(N+1))=-2*Vn*p_v;                            % deltaV (doppio prodotto)
end


pesiH((nx_tilde4*(N+1)+1):nu:(nx_tilde4*(N+1)+nu*N))    =   p_u1_dp;
pesiH((nx_tilde4*(N+1)+2):nu:(nx_tilde4*(N+1)+nu*N))    =   p_u2_dp;
pesiH((nx_tilde4*(N+1)+3):nu:(nx_tilde4*(N+1)+nu*N))    =   p_u3_dp;
pesiH((nx_tilde4*(N+1)+4):nu:(nx_tilde4*(N+1)+nu*N))    =   p_u4_dp;

pesiH((nx_tilde4*(N+1)+5):nu:(nx_tilde4*(N+1)+nu*N))    =   p_u1_dq;
pesiH((nx_tilde4*(N+1)+6):nu:(nx_tilde4*(N+1)+nu*N))    =   p_u2_dq;
pesiH((nx_tilde4*(N+1)+7):nu:(nx_tilde4*(N+1)+nu*N))    =   p_u3_dq;
pesiH((nx_tilde4*(N+1)+8):nu:(nx_tilde4*(N+1)+nu*N))    =   p_u4_dq;

pesiH((nx_tilde1+1):nx_tilde4:nx_tilde4*(N+1))          =    p_u1_p;
pesiH((nx_tilde1+2):nx_tilde4:nx_tilde4*(N+1))          =    p_u2_p;
pesiH((nx_tilde1+3):nx_tilde4:nx_tilde4*(N+1))          =    p_u3_p;
pesiH((nx_tilde1+4):nx_tilde4:nx_tilde4*(N+1))          =    p_u4_p;

pesiH((nx_tilde1+5):nx_tilde4:nx_tilde4*(N+1))          =    p_u1_q;
pesiH((nx_tilde1+6):nx_tilde4:nx_tilde4*(N+1))          =    p_u2_q;
pesiH((nx_tilde1+7):nx_tilde4:nx_tilde4*(N+1))          =    p_u3_q;
pesiH((nx_tilde1+8):nx_tilde4:nx_tilde4*(N+1))          =    p_u4_q;


H_=diag(pesiH);
H_=2*H_;

[Output, Jopt, exitflag, outputsolver] = cplexqp(H_, f_, Aineq_, bineq_, Aeq_, beq_);

%In order to interrupt simulation only if feasibility is compromised
if exitflag~=1 && exitflag~=5
    disp('Errore solutore')
    keyboard
end

% Saving of Output Variables of the MPC Controller

Umpc            =       reshape(Output(nx_tilde4*(N+1)+1:nx_tilde4*(N+1)+nu*N),nu,N);
Xmpc            =       reshape(Output(1:nx_tilde4*(N+1)),nx_tilde4,N+1);
v               =       Xmpc(2*n_nodes+1,2);
V               =       Xmpc(1:n_nodes,2);
Soc             =       Xmpc((nx_tilde3+1) : (nx_tilde4),2);
TotU            =       Xmpc((nx_tilde1+1) : nx_tilde2,2);
Pgen            =       Xmpc((nx_tilde2+1) : nx_tilde3,2);
val_h           =       Output(i_h);
val_hpsoft      =       Output(i_hpsoft);

end

 