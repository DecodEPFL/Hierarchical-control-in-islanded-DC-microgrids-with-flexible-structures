if  (Pred_PV_k) < 0        
    
    sourcenodes_sim = [sourcenodes pvnodes];      

    Vref_sim = Vref;

else
    
    sourcenodes_sim = [sourcenodes];
    
    Vref_sim = Vref(1:end-1);
    
end

    
    I_L_sim = blkdiag(zeros(n+p),I_L(p+1:end,p+1:end))*ones(nn,1);
    P_L_sim = blkdiag(zeros(n+p),P_L(p+1:end,p+1:end));
    Y_L_sim = blkdiag(zeros(n+p),Y_L(p+1:end,p+1:end));
    

G_s_sim = G_s;


v_sim            =  sdpvar(nn,1);
p_G_sim          =  sdpvar(nn,1);  
i_g_sim          =  sdpvar(nn,1);
p_gen_sim        =  sdpvar(n+p,1);

constr               = [ diag(v_sim)*Lap_t*v_sim + diag(v_sim)*I_L_sim + diag(i_g_sim)*G_s_sim*(R_t_mod)*G_s_sim'*i_g_sim +  P_L_sim*ones(nn,1) - G_s_sim*(p_gen_sim) ==  zeros(nn,1); ...
                         v_sim([sourcenodes_sim]) ==  Vref_sim; ...
                         v_sim([setdiff(1:nn,sourcenodes_sim)]) <= 1.2*v_nom([setdiff(1:nn,sourcenodes_sim)]); ...
                         v_sim([setdiff(1:nn,sourcenodes_sim)]) >= 0.8*v_nom([setdiff(1:nn,sourcenodes_sim)]); ...
                       % p_G_sim   ==  (diag(v_sim))*i_g_sim +  diag(i_g_sim)*G_s_sim*(R_t_mod)*G_s_sim'*i_g_sim; ...  
                         (diag(v_sim))*i_g_sim  ==  (diag(v_sim))*Lap*v_sim + (diag(v_sim))*I_L_sim  + (diag(v_sim))*Y_L_sim*v_sim + P_L_sim*ones(nn,1)]; 

if  (Pred_PV_k) == 0   
    constr           = [ constr;...
                         p_gen_sim(end) == -P_L(1,1);];
end



J  = norm(v_sim-v_opf,2);


assign(v_sim,v_opf);
 assign(p_G_sim,double(p_G_opt));
assign(i_g_sim,double(i_g_opt));

opt=sdpsettings('solver','fmincon','showprogress',0,'verbose',0,'usex0',1);
opt.fmincon.MaxFunEvals= 30000;
opt.fmincon.MaxIter= 30000;
opt.fmincon.TolCon = 1e-3;
sol=solvesdp(constr,J,opt);
ii=0;

if and(sol.problem==3,t>1)   
    keyboard
end

%%
while and(sol.problem==1,ii<10)
   disp('INFEASIBLE'); 
   
   ii= ii+1;

assign(v_sim,+1*v_opf);
opt=sdpsettings('solver','fmincon','showprogress',1,'verbose',1,'usex0',1);
opt.fmincon.MaxFunEvals= 30000;
opt.fmincon.MaxIter= 30000;
opt.fmincon.TolCon = 1e-3;
sol=solvesdp(constr,J,opt);

    
if ii ==10
    keyboard;
end

end
 
v_sim_opt             = double(v_sim);
p_G_sim_opt           = double(p_gen_sim);
