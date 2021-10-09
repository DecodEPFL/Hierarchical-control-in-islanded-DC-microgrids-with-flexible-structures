function [J,Pt,Qt,Y] = RealJacobian(V,angleV,w,RegG,RegL,netdata)

%% Network data

ni  =  netdata.ni;
nf  =  netdata.nf;
R   =  netdata.R;
L   =  netdata.L;

n_nodes    =  netdata.n_nodes;
nodi_gen   =  netdata.nodi_gen;
nodi_load  =  netdata.nodi_load;

%% Impedance matrix                                                         

X   =  w*L;  


Yreal  =  zeros(n_nodes,n_nodes);                                              
Yimag  =  zeros(n_nodes,n_nodes);                                              

dyrealdw  =  zeros(n_nodes,n_nodes);
dyimagdw  =  zeros(n_nodes,n_nodes);

% Definition of extra-diagonal terms

for k=1:length(ni)                                                        
    a=ni(k);
    b=nf(k);
    
    Yreal(a,b) =    - R(k) / ( R(k)^2 + X(k)^2 );                           
    Yreal(b,a) =    Yreal(a,b);                                            % the matrix is symmetric
    
    Yimag(a,b) =    + X(k) / ( R(k)^2 + X(k)^2 );                           
    Yimag(b,a) =    Yimag(a,b);                                            % the matrix is symmetric
    
    
    % Derivatives computation
    
    dyrealdw(a,b) =  R(k)*(2*w*L(k)^2)/ (R(k)^2 + X(k)^2)^(2);                                        
    dyrealdw(b,a) =  dyrealdw(a,b);
    
    dyimagdw(a,b) =  (L(k)*(R(k)^2 + X(k)^2) - L(k)*w*(2*w*L(k)*L(k)))/(R(k)^2 + X(k)^2)^2;
    dyimagdw(b,a) =  dyimagdw(a,b);
    
end

% Definition of diagonal terms
for k=1:n_nodes                                                            
    
     for j=1:n_nodes
        
         if(j~=k)
            Yreal(k,k) = Yreal(k,k) - Yreal(k,j);                          
            Yimag(k,k) = Yimag(k,k) - Yimag(k,j);
            
            dyrealdw(k,k) = dyrealdw(k,k) - dyrealdw(k,j);
            dyimagdw(k,k) = dyimagdw(k,k) - dyimagdw(k,j);
         end
        
     end
end

Y=Yreal + 1i*Yimag;                                                        

%% Trasmitted powers between nodes


Pt=zeros(1,n_nodes);                                                   
Qt=zeros(1,n_nodes);

for i=1:n_nodes
        
    for j=1:n_nodes
                        
        if(j~=i)
                        Pt(i)= Pt(i) + V(i)*V(j)*(  ( Yreal(i,j) )*cos(angleV(i)-angleV(j))   +   ( Yimag(i,j) )*sin(angleV(i)-angleV(j)) );
                        Qt(i)= Qt(i) + V(i)*V(j)*(  ( Yreal(i,j) )*sin(angleV(i)-angleV(j))   -   ( Yimag(i,j) )*cos(angleV(i)-angleV(j)) );
        end    
                            
    end
    
    Pt(i)= Pt(i)+Yreal(i,i)*V(i)*V(i);
    Qt(i)= Qt(i)-Yimag(i,i)*V(i)*V(i);
    
end

%% Derivatives of trasmitted powers

diffPt_V            =  zeros(n_nodes,n_nodes);
diffPt_angleV       =  zeros(n_nodes,n_nodes);

diffQt_V            =  zeros(n_nodes,n_nodes);
diffQt_angleV       =  zeros(n_nodes,n_nodes);

diffPt_w            =  zeros(1,n_nodes);
diffQt_w            =  zeros(1,n_nodes);

for i=1:n_nodes
        
        for j=1:n_nodes
            
            
            if (j~=i)
                
            diffPt_V(i,i)= diffPt_V(i,i) + V(j)*(  ( Yreal(i,j) )*cos(angleV(i)-angleV(j))   +   ( Yimag(i,j) )*sin(angleV(i)-angleV(j)) );
            diffPt_V(i,j)= V(i)*(  ( Yreal(i,j) )*cos(angleV(i)-angleV(j))   +   ( Yimag(i,j) )*sin(angleV(i)-angleV(j)) );
            
            
            diffQt_V(i,i)=diffQt_V(i,i) + V(j)*(  ( Yreal(i,j) )*sin(angleV(i)-angleV(j))   -   ( Yimag(i,j) )*cos(angleV(i)-angleV(j)) );
            diffQt_V(i,j)= V(i)*(  ( Yreal(i,j) )*sin(angleV(i)-angleV(j))   -   ( Yimag(i,j) )*cos(angleV(i)-angleV(j)) );
            
                        
            diffPt_angleV(i,i)= diffPt_angleV(i,i) + V(i)*V(j)*(  ( Yreal(i,j) )*(-sin(angleV(i)-angleV(j)))   +   ( Yimag(i,j) )*cos(angleV(i)-angleV(j)) );
            diffPt_angleV(i,j)= diffPt_angleV(i,j) + V(i)*V(j)*(  ( Yreal(i,j) )*sin(angleV(i)-angleV(j))      +   ( Yimag(i,j) )*(-cos(angleV(i)-angleV(j))));
                     
            
            diffQt_angleV(i,i)= diffQt_angleV(i,i) + V(i)*V(j)*(  ( Yreal(i,j) )*(cos(angleV(i)-angleV(j)))     +   ( Yimag(i,j) )*sin(angleV(i)-angleV(j)) );
            diffQt_angleV(i,j)= diffQt_angleV(i,j) + V(i)*V(j)*(  ( Yreal(i,j) )*(-cos(angleV(i)-angleV(j)))    -   ( Yimag(i,j) )*(sin(angleV(i)-angleV(j))));
            
            
            diffPt_w(1,i) = diffPt_w(1,i) + V(i)*V(j)*((dyrealdw(i,j))*cos(angleV(i)-angleV(j))   +   ( dyimagdw(i,j) )*sin(angleV(i)-angleV(j)) );
            
            diffQt_w(1,i) = diffQt_w(1,i) + V(i)*V(j)*((dyrealdw(i,j))*sin(angleV(i)-angleV(j))   -   ( dyimagdw(i,j) )*cos(angleV(i)-angleV(j)) );

            end
            
            
            
        end

        diffPt_V(i,i) = diffPt_V(i,i) +2*V(i)*Yreal(i,i);
        diffQt_V(i,i) = diffQt_V(i,i) +2*V(i)*(-Yimag(i,i));
        
        diffPt_w(1,i) = diffPt_w(1,i) + V(i)*V(i)*dyrealdw(i,i);
        diffQt_w(1,i) = diffQt_w(1,i) + V(i)*V(i)*(-dyimagdw(i,i));
        
end
                         
%% Definition of loads' powers

Psl     =  zeros(1,n_nodes);           
Qsl     =  zeros(1,n_nodes);
dPsldV  =  zeros(1,n_nodes);
dPsldw  =  zeros(1,n_nodes);
dQsldV  =  zeros(1,n_nodes);
dQsldw  =  zeros(1,n_nodes);

for i=1:n_nodes
    
        for j=1:length(nodi_load)
            
            if( i == nodi_load(j))
                [Psl(i),Qsl(i),dPsldV(i),dPsldw(i),dQsldV(i),dQsldw(i)]=regolazione_carico(V(i),w,RegL(j));
            end
            
        end
    
end

%% Definition of generators' powers

Psg    = zeros(1,n_nodes);
Qsg    = zeros(1,n_nodes);
dPsgdV = zeros(1,n_nodes);
dPsgdw = zeros(1,n_nodes);
dQsgdV = zeros(1,n_nodes);
dQsgdw = zeros(1,n_nodes);


for i=1:n_nodes
        for j=1:length(nodi_gen)
            if(i==nodi_gen(j))
               [Psg(i),Qsg(i),dPsgdV(i),dPsgdw(i),dQsgdV(i),dQsgdw(i)]=regolazione_generatore(V(i),w,RegG(j));
            end
        end
end

%% Jacobian computation

Jp_v      =  zeros(n_nodes,n_nodes);
Jp_angleV =  diffPt_angleV(1:n_nodes,2:n_nodes);
Jp_w      =  zeros(n_nodes,1);

Jq_v      =  zeros(n_nodes,n_nodes);
Jq_angleV =  diffQt_angleV(1:n_nodes,2:n_nodes);
Jq_w      =  zeros(n_nodes,1);


for i=1:n_nodes
    
            Jp_w(i,1) = diffPt_w(i)   - dPsldw(i)  - dPsgdw(i);
            Jq_w(i,1) = diffQt_w(i)   - dQsldw(i)  - dQsgdw(i);
            

     for j=1:(n_nodes)
         
         if(i==j)                                                          
            Jp_v(i,j) = diffPt_V(i,j) - dPsldV(i)  - dPsgdV(i);
            Jq_v(i,j) = diffQt_V(i,j) - dQsldV(i)  - dQsgdV(i);
         else
            Jp_v(i,j) = diffPt_V(i,j);
            Jq_v(i,j) = diffQt_V(i,j);
         end
      
     end
       
         
end

J= [Jp_v, Jp_angleV, Jp_w; Jq_v, Jq_angleV, Jq_w];

end
