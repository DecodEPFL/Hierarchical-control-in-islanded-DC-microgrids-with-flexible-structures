%% Inizalitation

close all

clear
clc

yalmip('clear') 
clear classes

Initialization;


%%
Umpc = zeros(n_b*2,1);

% S_load(30:50,2) = S_load(30:50,2) + 0.5*Srif*ones(21,1);
% S_load(120:150,2) = S_load(120:150,2) + 0.5*Srif*ones(31,1);
% S_load(70:95,2) = S_load(70:95,2) + 1i*0.5*Srif*ones(26,1);

for t = 1 : (Tsim)
    
    results = connected_PF(netdata,RegG,RegL);
   
     t
     V           =  abs(results.V);
     angleV      =  angle(results.V);
        
    [J,Pt,Qt,Y] = RealJacobian(V, angleV, netdata);
    Jinv  =  inv(J);
    
    [Jloss, Ploss] = LossesJacobian(V, angleV, netdata);
    
    P_batt(1) = Pt(5)-real(S_gen(t,1))/Srif;
    P_batt(2) = Pt(6);
    Q_batt(1) = Qt(5);
    Q_batt(2) = Qt(6);
    Ps        = Pt(1);
    Qs        = Qt(1);
    SOC_pq(1) = SOC(1);
    SOC_pq(2) = SOC(2); 
    SOC_s     = SOC(3);
    
    
    x = [V(2:end); angleV(2:end); Ps; Qs; P_batt(1); P_batt(2); Q_batt(1); Q_batt(2); SOC_pq(1); SOC_pq(2); SOC_s];
    
    [Umpc, Jopt, h, SOC] = MPC_control(x, Jinv, Jloss, Ploss, var_D(:,(t):(t+N-1)), netdata, MPCpar, Umpc,t);
  
  % Aggiornamento Generatori
  
    RegG(1).P   = Umpc(1,1)  +  P_batt(1) + real(S_gen(t+1,1))/Srif;
    RegG(1).Q   = Umpc(1+length(nodi_ctrl_gen),1)  +  Q_batt(1);
    RegG(2).P   = Umpc(2,1)  +  P_batt(2);
    RegG(2).Q   = Umpc(2+length(nodi_ctrl_gen),1)  +  Q_batt(2);

  % Aggiornamento Carichi
 
    RegL(1).P  = -real(S_load(t+1,1))/Srif;
    RegL(1).Q  = -imag(S_load(t+1,1))/Srif;

    RegL(2).P  = -real(S_load(t+1,2))/Srif;
    RegL(2).Q  = -imag(S_load(t+1,2))/Srif;

  % Salvataggio Dati    
    datalog.P_batt_1(t+1) = Umpc(1,1)  +  P_batt(1);
    datalog.P_node_5(t+1) = Umpc(1,1)  +  P_batt(1) + real(S_gen(t+1,1))/Srif;
    datalog.P_batt_2(t+1) = Umpc(2,1)  +  P_batt(2);    
    datalog.Q_batt_1(t+1) = Umpc(1+length(nodi_ctrl_gen),1)  +  Q_batt(1);
    datalog.Q_batt_2(t+1) = Umpc(2+length(nodi_ctrl_gen),1)  +  Q_batt(2);  
    
    datalog.Umpc(:,:,t)   =  Umpc;
    datalog.Jopt(t)       =  Jopt;
    datalog.h(:,:,t)      =  h;
    
    datalog.x(:,t)        =  x;
    datalog.SOC(:,t+1)    =  SOC;
    datalog.P_loss(:,t)   =  Ploss; 
end

t=Tsim;

    results = connected_PF(netdata,RegG,RegL);

    V           =  abs(results.V);
    angleV      =  angle(results.V);

    
   x = [V(2:end); angleV(2:end)];                                              %Definisco lo stato attuale
    datalog.x(1:size(x,1),t+1)        =  x;

 %%
 
 figure1 = figure(1);
 plot(MPCpar.step*60:MPCpar.step*60:h_islanded*60, datalog.x(1:n_nodes-1,1:Tsim),'LineWidth',2 );
 grid on
 ylim([0.9 1.1])
 title('Voltages [p.u.]')
 xlabel('Time [min]');
 
  
figure2 = figure(2);
stairs(MPCpar.step*60:MPCpar.step*60:h_islanded*60, datalog.P_batt_1(1,1:Tsim),'LineWidth',2 );
hold on
stairs(MPCpar.step*60:MPCpar.step*60:h_islanded*60, datalog.P_batt_2(1,1:Tsim),'LineWidth',2 );
hold on
plot(MPCpar.step*60:MPCpar.step*60:h_islanded*60, datalog.x(2*(n_nodes-1)+1,1:Tsim),'LineWidth',2 );
grid on
legend('Batt 1','Batt 2','Slack');
title('Generation active powers')
xlabel('Time [min]');

figure3 = figure(3);
stairs(MPCpar.step*60:MPCpar.step*60:h_islanded*60, datalog.Q_batt_1(1,1:Tsim),'LineWidth',2 );
hold on
stairs(MPCpar.step*60:MPCpar.step*60:h_islanded*60, datalog.Q_batt_2(1,1:Tsim),'LineWidth',2 );
hold on
plot(MPCpar.step*60:MPCpar.step*60:h_islanded*60, datalog.x(2*(n_nodes-1)+2,1:Tsim),'LineWidth',2 );
grid on
legend('Batt 1','Batt 2','Slack');
title('Generation reactive powers')
xlabel('Time [min]');

figure4 = figure(4);
plot(MPCpar.step*60:MPCpar.step*60:h_islanded*60, datalog.SOC(:,1:Tsim),'LineWidth',2 );
grid on
title('SOCs')
legend('Batt 1','Batt 2','Slack');
xlabel('Time [min]');

figure5 = figure(5);
plot(MPCpar.step*60:MPCpar.step*60:h_islanded*60, real(S_load(1:Tsim,1))/Srif,'LineWidth',2 );
hold on
plot(MPCpar.step*60:MPCpar.step*60:h_islanded*60, real(S_load(1:Tsim,2))/Srif,'LineWidth',2 );
hold on
plot(MPCpar.step*60:MPCpar.step*60:h_islanded*60, real(S_gen(1:Tsim,1))/Srif,'LineWidth',2 );
grid on
title('Non dispatchable active powers')
legend('Load 1','Load 2','PV');
xlabel('Time [min]');


figure6 = figure(6);
plot(MPCpar.step*60:MPCpar.step*60:h_islanded*60, imag(S_load(1:Tsim,1))/Srif,'LineWidth',2 );
hold on
plot(MPCpar.step*60:MPCpar.step*60:h_islanded*60, imag(S_load(1:Tsim,2))/Srif,'LineWidth',2 );
grid on
title('Non dispatchable reactive powers')
legend('Load 1','Load 2','PV');
xlabel('Time [min]');

figure7 = figure(7);
plot(MPCpar.step*60:MPCpar.step*60:h_islanded*60, +real(S_gen(1:Tsim,1)')/Srif + datalog.P_batt_1(1,1:Tsim) + datalog.P_batt_2(1,1:Tsim) + datalog.x(2*(n_nodes-1)+1,1:Tsim),'LineWidth',2 );
hold on
plot(MPCpar.step*60:MPCpar.step*60:h_islanded*60, real(S_load(1:Tsim,1))/Srif + real(S_load(1:Tsim,2))/Srif ,'LineWidth',2 );
grid on
legend('Total generation','Absorption');
title('Active powers')
xlabel('Time [min]');
% 
figure8 = figure(8);
plot(MPCpar.step*60:MPCpar.step*60:h_islanded*60, sum(datalog.P_loss(:,1:Tsim),1),'LineWidth',2 );
legend('No Plosses modelling','Ploss minimization');
grid
title('Total losses')
xlabel('Time [min]');
% 
% 
% % % 
% figure
% plot(MPCpar.step*60:MPCpar.step*60:h_islanded*60, datalog.Q_batt_1(1,1:Tsim) + datalog.Q_batt_2(1,1:Tsim) + datalog.x(2*(n_nodes-1)+2,1:Tsim),'LineWidth',2 );
% hold on
% plot(MPCpar.step*60:MPCpar.step*60:h_islanded*60, imag(S_load(1:Tsim,1))/Srif + imag(S_load(1:Tsim,2))/Srif ,'LineWidth',2 );
% grid on
% legend('Total generation','Absorption');
% title('Reactive powers')
% xlabel('Time [min]');
