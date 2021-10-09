function [J,Pt,Qt,Y] = RealJacobian(V,angleV,netdata)

%% Network data

ni  =  netdata.ni;
nf  =  netdata.nf;
R   =  netdata.R;
L   =  netdata.L;

n_nodes    =  netdata.n_nodes;
% nodi_load  =  netdata.nodi_load;

w=2*50*pi;

%% Impedance matrix                                                         

X   =  w*L;  

Yreal  =  zeros(n_nodes,n_nodes);                                              
Yimag  =  zeros(n_nodes,n_nodes);                                             

% Definition of extra-diagonal terms

for k=1:length(ni)                                                        
    a=ni(k);
    b=nf(k);
    
    Yreal(a,b) =    - R(k) / ( R(k)^2 + X(k)^2 );                           
    Yreal(b,a) =    Yreal(a,b);                                            % the matrix is symmetric
    
    Yimag(a,b) =    + X(k) / ( R(k)^2 + X(k)^2 );                           
    Yimag(b,a) =    Yimag(a,b);                                            % the matrix is symmetric
  
end

% Definition of diagonal terms
for k=1:n_nodes                                                               
     for j=1:n_nodes  
         if(j~=k)
            Yreal(k,k) = Yreal(k,k) - Yreal(k,j);                          
            Yimag(k,k) = Yimag(k,k) - Yimag(k,j);
         end
        
     end
end

for i=1:6

    Yimag(i,i) = Yimag(i,i);
end

Y = Yreal + 1i*Yimag;                                                        

%% Trasmitted powers between nodes


Pt=zeros(1,n_nodes);                                                   
Qt=zeros(1,n_nodes);

% for i=1:n_nodes
%         
%     for j=1:n_nodes
%                         
%         if(j~=i)
%                         Pt(i)= Pt(i) + V(i)*V(j)*(  ( Yreal(i,j) )*cos(angleV(i)-angleV(j))   +   ( Yimag(i,j) )*sin(angleV(i)-angleV(j)) );
%                         Qt(i)= Qt(i) + V(i)*V(j)*(  ( Yreal(i,j) )*sin(angleV(i)-angleV(j))   -   ( Yimag(i,j) )*cos(angleV(i)-angleV(j)) );
%         end    
%                             
%     end
%     
%     Pt(i)= Pt(i)+Yreal(i,i)*V(i)*V(i);
%     Qt(i)= Qt(i)-Yimag(i,i)*V(i)*V(i);
%     
% end


for i=1:n_nodes
        
    for j=1:n_nodes
                        
                        Pt(i)= Pt(i) +  V(i)*V(j)*(abs(Y(i,j))*cos(  angleV(i)-angleV(j) - angle(Y(i,j))  ) );
                        Qt(i)= Qt(i) +  V(i)*V(j)*abs(Y(i,j))*sin(  angleV(i)-angleV(j) - angle(Y(i,j))  );
   
                            
    end
    

end

%% Derivatives of trasmitted powers

diffPt_V            =  zeros(n_nodes,n_nodes);
diffPt_angleV       =  zeros(n_nodes,n_nodes);

diffQt_V            =  zeros(n_nodes,n_nodes);
diffQt_angleV       =  zeros(n_nodes,n_nodes);

for i=1:n_nodes 
    
        for j=1:n_nodes
            
          if (j~=i)
                
            diffPt_V(i,i)= diffPt_V(i,i) + V(j)*(abs(Y(i,j)))*cos(angleV(i)-angleV(j) - angle(Y(i,j)));
            diffPt_V(i,j)= V(i)*((abs(Y(i,j)))*cos(angleV(i)-angleV(j)-angle(Y(i,j))));
            
             
            diffQt_V(i,i) = diffQt_V(i,i) + V(j)*(abs(Y(i,j)))*sin(angleV(i)-angleV(j) - angle(Y(i,j)));
            diffQt_V(i,j) = V(i)*(abs(Y(i,j)))*sin(angleV(i)-angleV(j)-angle(Y(i,j)));
            
                        
            diffPt_angleV(i,i)= diffPt_angleV(i,i) - V(i)*V(j)*(abs(Y(i,j)))*sin(angleV(i)-angleV(j)-angle(Y(i,j)));
            diffPt_angleV(i,j)=  V(i)*V(j)*(abs(Y(i,j)))*sin(angleV(i)-angleV(j)-angle(Y(i,j)));
                     
            
            diffQt_angleV(i,i)= diffQt_angleV(i,i) + V(i)*V(j)*(abs(Y(i,j)))*cos(angleV(i)-angleV(j)-angle(Y(i,j)));
            diffQt_angleV(i,j)= - V(i)*V(j)*(abs(Y(i,j)))*cos(angleV(i)-angleV(j)-angle(Y(i,j)));            
            end      
        end

        diffPt_V(i,i) = diffPt_V(i,i) + 2*V(i)*abs(Y(i,i))*cos(angle(Y(i,i)));
        diffQt_V(i,i) = diffQt_V(i,i) - 2*V(i)*abs(Y(i,i))*sin(angle(Y(i,i)));         
end
                         
%% Definition of loads' powers

Psl     =  zeros(1,n_nodes);           
Qsl     =  zeros(1,n_nodes);
dPsldV  =  zeros(1,n_nodes);
dQsldV  =  zeros(1,n_nodes);
% 
% for i=1:n_nodes
%     
%         for j=1:length(nodi_load)
%             
%             if( i == nodi_load(j))
%                 [Psl(i),Qsl(i),dPsldV(i),dPsldw(i),dQsldV(i),dQsldw(i)]=regolazione_carico(V(i),w,RegL(j));
%             end
%             
%         end
%     
% end

%% Definition of generators' powers

Psg    = zeros(1,n_nodes);
Qsg    = zeros(1,n_nodes);
dPsgdV = zeros(1,n_nodes);
dQsgdV = zeros(1,n_nodes);

% 
% for i=1:n_nodes
%         for j=1:length(nodi_gen)
%             if(i==nodi_gen(j))
%                [Psg(i),Qsg(i),dPsgdV(i),dPsgdw(i),dQsgdV(i),dQsgdw(i)]=regolazione_generatore(V(i),w,RegG(j));
%             end
%         end
% end

%% Jacobian computation

Jp_v      =  diffPt_V(2:end,2:end);
Jp_angleV =  diffPt_angleV(2:end,2:end);

Jq_v      =  diffQt_V(2:end,2:end);
Jq_angleV =  diffQt_angleV(2:n_nodes,2:n_nodes);


J = [Jp_v, Jp_angleV; Jq_v, Jq_angleV];

end
