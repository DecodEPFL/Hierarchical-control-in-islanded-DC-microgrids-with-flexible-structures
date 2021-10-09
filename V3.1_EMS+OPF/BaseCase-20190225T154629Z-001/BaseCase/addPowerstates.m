function [A_tilde3, B_tilde3, M_tilde3, nx_tilde3, Dprof_tilde3] = addPowerstates(A_tilde2, B_tilde2, M_tilde2, RegG, N, nodi_ctrl_gen, n_nodes, x, nx_tilde2, nx_tilde1, Dprof_tilde2, Vref)

ncg=length(nodi_ctrl_gen);                                               %Number of power states

m_v=zeros(1,ncg);
m_w=zeros(1,ncg);

for j = 1:ncg
    [~,~,dPdV,~,~,dQdw]=regolazione_generatore(x(nodi_ctrl_gen(j)),x(2*n_nodes),RegG(j));
    
     m_v(j)=dPdV;
     m_w(j)=dQdw;       
end

%%% I add to the actual states the estiamated generated active and reactive
%%% powers

nx_tilde3=nx_tilde2 +2*ncg;
%% Creation A_tilde3

%%% The matrix A_tilde3 is created by plugging different blocks

A_tilde3=zeros(nx_tilde3,nx_tilde3);

A_tilde3(1:nx_tilde2,1:nx_tilde2) = A_tilde2;

A_tilde3((nx_tilde2+1):nx_tilde3,nx_tilde1+1:nx_tilde2) = eye(2*ncg);

for j=1:ncg
   A_tilde3((nx_tilde2)+j,nodi_ctrl_gen(j))   =   m_v(j);
   A_tilde3((nx_tilde2)+j+ncg,2*n_nodes)      =   m_w(j);
end



%% Creation of B_tilde3

Bcap1=eye(2*ncg);
Bcap2=zeros(2*ncg,2*ncg);

for j=1:ncg
   Bcap2(j,1:2*ncg)=m_v(j)*B_tilde2(nodi_ctrl_gen(j),1:2*ncg);
   Bcap2(j+ncg,1:2*ncg)=m_w(j)*B_tilde2(2*n_nodes,1:2*ncg);
end

Bcap=Bcap1+Bcap2;

B_tilde3=[B_tilde2;Bcap];

%% Creation of M_tilde3


Mcap1=zeros(2*ncg,2*n_nodes);

    for j=1:ncg
        Mcap1(j,1:end) = m_v(j)*M_tilde2(nodi_ctrl_gen(j),1:2*n_nodes);
        Mcap1(j+ncg,1:end) = m_w(j)*M_tilde2(2*n_nodes,1:2*n_nodes);
    end

Mcap2=zeros(2*ncg,2);

    for j=1:ncg
        Mcap2(j,2)=-m_v(j);
        Mcap2(j+ncg,1)=-m_w(j);
    end

M_tilde3=[M_tilde2, zeros(nx_tilde2,1);Mcap1, Mcap2];

Dprof_tilde3=[Dprof_tilde2; ones(1,size(Dprof_tilde2,2))*Vref];


end