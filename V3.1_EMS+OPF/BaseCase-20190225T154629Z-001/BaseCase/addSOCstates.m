function [A_tilde4, B_tilde4, M_tilde4, nx_tilde4, Dprof_tilde4] = addSOCstates(A_tilde3, B_tilde3, M_tilde3, gen_batt, nodi_ctrl_gen, nx_tilde3,k, Dprof_tilde3)

nb=length(gen_batt);
ncg=length(nodi_ctrl_gen);

nx_tilde4=nx_tilde3+nb;

%% Creation A_tilde4

Acap1=zeros(nb,nx_tilde3);

for b=1:nb
    
    for j=1:ncg
    
        if(nodi_ctrl_gen(j)==gen_batt(b))
        
            Acap1(b,:)=-k(b)*A_tilde3(j+nx_tilde3-2*ncg,:);
        
        end
    end
end

A_tilde4=[A_tilde3, zeros(nx_tilde3,nb); Acap1, eye(nb,nb)];



%% Creation B_tilde4
Bcap1=zeros(nb,2*ncg);

for b=1:nb
    
    for j=1:ncg
    
        if(nodi_ctrl_gen(j)==gen_batt(b))
        
            Bcap1(b,:)=-k(b)*B_tilde3(j+nx_tilde3-2*ncg,:);
        
        end
    end
end

B_tilde4=[B_tilde3  ;   Bcap1];

%% Creation M_tilde4

Mcap1=zeros(nb,size(M_tilde3,2));

for b=1:nb
    
    for j=1:ncg
    
        if(nodi_ctrl_gen(j)==gen_batt(b))
        
            Mcap1(b,:)=-k(b)*M_tilde3(j+nx_tilde3-2*ncg,:);
        
        end
    end
end

M_tilde4=[M_tilde3  ;   Mcap1];
%% Dprof

Dprof_tilde4 = Dprof_tilde3;
end