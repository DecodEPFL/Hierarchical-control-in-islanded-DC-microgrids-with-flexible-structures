clear
clc
close all

load profiles

dati_rete;

for t=1:size(Pload,1);
    
    SL1=Pload(t,1);
    SL2=Pload(t,2);
    SL3=Pload(t,3);
    SL4=Pload(t,4);
    SL5=Pload(t,5);
    
    PG1=Pgen(t,1);
    PG2=Pgen(t,2);
    PG3=Pgen(t,3);
    
    generatori_e_carichi;
    [results,details]=connected_PF(netdata,RegG,RegL);
    
    V(t,:)=abs(results.V);
    I(t,:)=abs(results.I);
    S(t,:)=details.Sgen;
%     f(t)=results.w/(2*pi);
    

end

figure(2)
plot(V(:,2:end),'linewidth',1.5)
ylim([360,440]);
grid on

figure(3)
plot(real(S),'linewidth',1.5)
grid on


