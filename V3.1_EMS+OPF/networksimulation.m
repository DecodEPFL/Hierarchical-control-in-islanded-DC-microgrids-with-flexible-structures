% new filter parameters 

pvnodes_nongen = [];

if (Pred_PV_k) < 0       
        ng_mod = n + p;
        nl_mod = m;
        I_L_sim = blkdiag(zeros(ng_mod),I_L(2:end,2:end))*ones(nn,1);
        P_L_sim = blkdiag(zeros(ng_mod),P_L(2:end,2:end));
        Y_L_sim = blkdiag(zeros(ng_mod),Y_L(2:end,2:end));
 
        alpha_mod = alpha;
        beta_mod  = beta;
        gamma_mod = gamma;
        C_t_mod = C_t;
        C_loads_mod = C_loads;
        Vref_mod = Vref;   
else
        ng_mod = n;
        nl_mod = m + p; 
        I_L_sim = blkdiag(zeros(ng_mod),I_L)*ones(nn,1);
        P_L_sim = blkdiag(zeros(ng_mod),P_L);
        Y_L_sim = blkdiag(zeros(ng_mod),Y_L);
        
        alpha_mod = alpha(1:end-1,1:end-1);
        beta_mod  = beta(1:end-1,1:end-1);
        gamma_mod = gamma(1:end-1,1:end-1);
        C_t_mod = C_t(1:end-1,1:end-1);
        C_loads_mod = blkdiag(C_t(end,end),C_loads);
        Vref_mod = Vref(1:end-1);  
        
        pvnodes_nongen = pvnodes;

end

deactivated_nodes = [];
deactivated_lines = [];

nn_mod = nn;
mm_mod = mm;


for kk = 1 : sn
    if delta_g_k(kk) == 0
        
        ng_mod  = ng_mod-1;
        nn_mod  = nn_mod-1;
        mm_mod  = mm_mod - length(find(B(switchnodes(kk),:)));
        
        I_L_sim = I_L_sim(2:end);
        P_L_sim = P_L_sim(2:end,2:end);
        Y_L_sim = Y_L_sim(2:end,2:end);
        
        alpha_mod  =   alpha_mod([setdiff(1:size(alpha_mod,1),switchnodes(kk))],[setdiff(1:size(alpha_mod,2),switchnodes(kk))]);
        beta_mod   =   beta_mod([setdiff(1:size(beta_mod,1),switchnodes(kk))],[setdiff(1:size(beta_mod,2),switchnodes(kk))]);
        gamma_mod  =   gamma_mod([setdiff(1:size(gamma_mod,1),switchnodes(kk))],[setdiff(1:size(gamma_mod,2),switchnodes(kk))]);
        C_t_mod    =   C_t([setdiff(1:size(C_t_mod,1),switchnodes(kk))],[setdiff(1:size(C_t_mod,1),switchnodes(kk))]);
        Vref_mod   =   Vref([setdiff(1:ng_mod+1,switchnodes(kk))]);
        
        deactivated_nodes = [deactivated_nodes, switchnodes(kk)];
        deactivated_lines = [deactivated_lines, find(B(switchnodes(kk),:))];
    end
end

%inval=v_nom(1)*ones(nn_mod + 2*(ng_mod) + mm_mod,1);

B_mod =  B(([setdiff(1:size(B,1),deactivated_nodes)]),([setdiff(1:size(B,2),deactivated_lines)]));
R_mod =  R(([setdiff(1:size(R,1),deactivated_lines)]),([setdiff(1:size(R,2),deactivated_lines)]));
L_mod =  L(([setdiff(1:size(L,1),deactivated_lines)]),([setdiff(1:size(L,2),deactivated_lines)]));

Ct = blkdiag(C_t_mod,  C_loads_mod);

if t>0
    inval = yout(end,[setdiff(1:nn,deactivated_nodes),setdiff(nn+1:nn+n+p, nn + [pvnodes_nongen, deactivated_nodes]), setdiff(nn+n+p+1 : nn +2*(n+p),nn + n + p +[pvnodes_nongen, deactivated_nodes] ),setdiff(nn +2*(n+p) + 1 : nn +2*(n+p) + mm, nn +2*(n+p) + deactivated_lines ) ]);
else
    inval=v_nom(1)*ones(nn + 2*(n) + mm,1);
end
PL = -(Ct^-1)*P_L_sim;

X=[-Ct^(-1)*Y_L_sim,[C_t_mod^(-1);zeros(nl_mod,ng_mod)],zeros(ng_mod+nl_mod,ng_mod),-Ct^(-1)*B_mod; [alpha_mod,zeros(ng_mod,nl_mod)],beta_mod,gamma_mod,zeros(ng_mod,mm_mod);[-eye(ng_mod),zeros(ng_mod,nl_mod)],zeros(ng_mod),zeros(ng_mod),zeros(ng_mod,mm_mod); L_mod^(-1)*B_mod',zeros(mm_mod,ng_mod),zeros(mm_mod,ng_mod),-L_mod^(-1)*R_mod];
I=[-Ct^(-1)*I_L_sim;zeros(ng_mod,1);Vref_mod;zeros(mm_mod,1)];

net_par.PL      = PL;
net_par.I       = I;
net_par.X       = X;
net_par.N       = nn_mod;

[time,yout_mod]=ode45(@(time,yout_mod) ZIPloads_par2(time,yout_mod,net_par),tr,inval);


yout = zeros(length(time), nn + 2*(n+p) + mm, 1);

yout(:,[setdiff(1:nn,deactivated_nodes)])                 =    yout_mod(:, 1:nn_mod);

yout(:, [setdiff(nn+1:nn+n+p, nn + [pvnodes_nongen, deactivated_nodes]) ])                          =    yout_mod(:, nn_mod+1:nn_mod + ng_mod);

yout(:, [setdiff(nn+n+p+1 : nn +2*(n+p),nn + n + p +[pvnodes_nongen, deactivated_nodes] )])         =    yout_mod(:, nn_mod+ng_mod+1:nn_mod + 2*ng_mod);

yout(:, [setdiff(nn +2*(n+p) + 1 : nn +2*(n+p) + mm, nn +2*(n+p) + deactivated_lines )])            =    yout_mod(:, nn_mod + 2*ng_mod+1:nn_mod + 2*ng_mod + mm_mod);
