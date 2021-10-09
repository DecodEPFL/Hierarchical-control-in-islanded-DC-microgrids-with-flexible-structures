yalmip('clear')


delta_nodes(switchnodes) = delta_g_k;

for kk = 1 : sn
    if delta_g_k(kk) == 0
        lines_g = find(B(switchnodes(kk),:));
        delta_lines(lines_g) = 0;
    else
        lines_g = find(B(switchnodes(kk),:));
        delta_lines(lines_g) = 1;
    end
end

I_L_opf = blkdiag(zeros(n+p),I_L(p+1:end,p+1:end))*ones(nn,1);
P_L_opf = blkdiag(zeros(n+p),P_L(p+1:end,p+1:end));
Y_L_opf = blkdiag(zeros(n+p),Y_L(p+1:end,p+1:end));

G_s  = zeros(nn,n+p);

G_s([sourcenodes pvnodes],1:n+p) = eye(n+p);

p_G  =  [P_ref_MPC; -P_L(1,1)];

%The Laplacian considers also if a line is disconnected or not
Lap    =  B*diag(delta_lines)*inv(R)*B';
Lap_t  =  Lap + Y_L_opf;

% If the generator is not connected the filter resistance is set to zero
R_t_mod                          =    R_t;
R_t_mod(switchnodes,switchnodes) =    diag(delta_g_k)*R_t(switchnodes,switchnodes);
  
Delta_P_slack    =  sdpvar(n+p,1); 
% 
if  (Pred_PV_k) < 0       
        Delta_P_slack_max_mod = [Delta_P_slack_max;  -Pred_PV_k];
        Delta_P_slack_min_mod = [Delta_P_slack_min;    P_L(1,1)];
else
        Delta_P_slack_max_mod = [Delta_P_slack_max;        1e-2];
        Delta_P_slack_min_mod = [Delta_P_slack_min;       -1e-2];
end


v_opt            =  sdpvar(nn,1);
p_G_opt          =  sdpvar(nn,1);  
i_g_opt          =  sdpvar(nn,1);

constr               = [ diag(v_opt)*Lap_t*v_opt + diag(v_opt)*I_L_opf + diag(i_g_opt)*G_s*(R_t_mod)*G_s'*i_g_opt +  P_L_opf*ones(nn,1) - G_s*(p_G + Delta_P_slack) ==  zeros(nn,1); ...
                         v_opt >=  0.9*v_nom.*delta_nodes ; ...
                         v_opt <=  1.1*v_nom.*delta_nodes ; ...
                         Delta_P_slack(:) <= Delta_P_slack_max_mod+1e-2;...
                         Delta_P_slack(:) >= Delta_P_slack_min_mod-1e-2;...                        
                         p_G_opt   ==  (diag(v_opt))*i_g_opt +  diag(i_g_opt)*G_s*(R_t_mod)*G_s'*i_g_opt; ...  
                        (diag(v_opt))*i_g_opt  ==  (diag(v_opt))*Lap*v_opt + (diag(v_opt))*I_L_opf  + (diag(v_opt))*Y_L_opf*v_opt + P_L_opf*ones(nn,1) ;...
                      ]; 


J  = norm(Delta_P_slack,2);


assign(v_opt,v_opf);
opt=sdpsettings('solver','fmincon','showprogress',0,'verbose',0,'usex0',1);
opt.fmincon.MaxFunEvals= 30000;
opt.fmincon.MaxIter= 30000;
%opt.fmincon.TolCon = 1e-4;

sol=solvesdp(constr,J,opt);

% if and(sol.problem==3,t>1)   
 %     keyboard
% end

ii=0;
%%
while and(or(sol.problem==1,sol.problem==3),ii<10)
   disp('OPF repeat');
   sol.problem
   
   ii= ii+1;

assign(v_opt,+1*v_nom);
opt=sdpsettings('solver','fmincon','showprogress',0,'verbose',0,'usex0',1);
opt.fmincon.MaxFunEvals= 30000;
opt.fmincon.MaxIter= 30000;
sol=solvesdp(constr,J,opt);
    
if ii ==10
    keyboard;
end

end
 
v_opf             = double(v_opt);
Delta_P_slack_opf = double(Delta_P_slack);


