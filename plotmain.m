%%
close all

load('data_HMPC_LoadMeasure.mat')

%%

P_real = zeros(p+m,Tsim);
vopf = datalog.OPF.v_opf;

for ii= 1:p
    P_real(ii,:) =  realProfiles.PV(ii).YL(:)'.*vopf(pvnodes(ii),:).^2 + realProfiles.PV(ii).I(:)'.*vopf(pvnodes(ii),:) + realProfiles.PV(ii).P(:)';    
end

for jj = 1:m
    P_real(p+jj,:) = realProfiles.L(jj).YL(:)'.*vopf(loadnodes(jj),:).^2 + realProfiles.L(jj).I(:)'.*vopf(loadnodes(jj),:) + realProfiles.L(jj).P(:)';     
end


%%

figure
plot((t_start:t_start+Tsim-1)/60,realProfiles.L(1).I, 'linewidth',2)
hold on
plot((t_start:t_start+Tsim-1)/60,realProfiles.L(4).I, 'linewidth',2)
hold on
plot((t_start:t_start+Tsim-1)/60,realProfiles.L(7).I, 'linewidth',2)
grid
xlim([0 24])
set(gca,'fontsize',18, 'FontName', 'Times New Roman')
ylabel('[A]');
ylim([0,0.065])
xlabel('Time [h]')
set(gca,'XTick',[0:4:24]);
legend({'L_A', 'L_B', 'L_C'},'Orientation','horizontal','location','north');
box on
%title('Loads constant current')

figure
plot((t_start:t_start+Tsim-1)/60,realProfiles.L(1).P, 'linewidth',2)
hold on
plot((t_start:t_start+Tsim-1)/60,realProfiles.L(4).P, 'linewidth',2)
hold on
plot((t_start:t_start+Tsim-1)/60,realProfiles.L(7).P, 'linewidth',2)
grid on 
set(gca,'fontsize',18, 'FontName', 'Times New Roman')
xlim([0 24])
ylabel('[kW]');
xlabel('Time [h]')
ylim([0.3,2])
legend({'L_A', 'L_B', 'L_C'},'Orientation','horizontal','location','north');
set(gca,'XTick',[0:4:24]);
box on
%title('Loads costant power')

%%

figure
%plot((t_start:t_start+Tsim-1)/60,P_nom(2,t_start : (t_start+Tsim)-1), 'linewidth',1)
hold on
plot((t_start:t_start+Tsim-1)/60,P_real(2, t_start : (t_start+Tsim)-1),'-r', 'linewidth',1.5)
hold on
plot((t_start:15:t_start+Tsim-1)/60,PL_nom_s(2,1 : 96), '--k', 'linewidth',2)
%title('Load power, type A')
grid on 
%legend({'Real','Forecast'},'Orientation','horizontal','location','north');
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
ylim([7,13])
xlim([0 24])
xlabel('Time [h]')
ylabel('[kW]');
set(gca,'XTick',[0:4:24]);
box on


figure
%plot((t_start:t_start+Tsim-1)/60,P_nom(5,t_start : (t_start+Tsim)-1), 'linewidth',1)
hold on
plot((t_start:t_start+Tsim-1)/60,P_real(5, t_start : (t_start+Tsim)-1),'-r', 'linewidth',1.5)
hold on
plot((t_start:15:t_start+Tsim-1)/60,PL_nom_s(5,1 : 96), '--k', 'linewidth',2)
%title('Load power, type B')
grid on 
%legend({'Real','Forecast'},'Orientation','horizontal','location','north');
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
xlim([0 24])
xlabel('Time [h]')
ylim([6.5,11])
ylabel('[kW]');
set(gca,'XTick',[0:4:24]);
box on

figure
%plot((t_start:t_start+Tsim-1)/60,P_nom(10,t_start : (t_start+Tsim)-1), 'linewidth',1)
hold on
plot((t_start:t_start+Tsim-1)/60,P_real(10, t_start : (t_start+Tsim)-1), '-r', 'linewidth',1.5)
hold on
plot((t_start:15:t_start+Tsim-1)/60,PL_nom_s(10,1 : 96), '--k', 'linewidth',2)
%title('Load power, type C')
grid on 
%legend({'Real','Forecast'},'Orientation','horizontal','location','north');
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
xlim([0 24])
ylim([8.3 16.3])
xlabel('Time [h]')
ylabel('[kW]');
set(gca,'XTick',[0:4:24]);
box on

%%
figure
plot((t_start:t_start+Tsim-1)/60, datalog.MPC.P_ref(end, t_start : (t_start+Tsim)-1)  , '--g', 'linewidth',2.5)
hold on
plot((t_start:t_start+Tsim-1)/60,datalog.OPF.Pg(pvnodes,1:end),'-r','linewidth',1.5)
hold on
plot((t_start:15:t_start+Tsim-1)/60,-P_nom(1,t_start :15: (t_start+Tsim)-1), ':','linewidth',1.5)
hold on
plot((t_start*MPCpar.t_s/60: MPCpar.t_s/60 :(t_start+Tsim-1)/60 ), -PL_nom_s(1, (t_start-1)*MPCpar.t_s + 1 : Tsim/MPCpar.t_s ), '--k', 'linewidth',1.5)
%legend({'Real','Nominal','Forecast'});
%title('PV system')
grid on 
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
xlim([0 24])
ylabel('[kW]');
set(gca,'XTick',[0:4:24]);
xlabel('Time [h]')
box on

%%
figure
%plot((t_start:t_start+Tsim-1)/60,100*datalog.SOC(1,1:end), ':','linewidth',2.5)
hold on
plot((t_start:t_start+Tsim-1)/60,datalog.OPFsim.p_G_sim(batt_nodes(1),1:end),'-r','linewidth',2.5)
hold on
plot((t_start:t_start+Tsim-1)/60,datalog.MPC.P_ref(batt_nodes(1),1:end), '--b','linewidth',2.5)
%title('Batt 1')
grid on 
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
%legend({'Real power','Reference power'},'Orientation','horizontal','location','north');
xlim([0 24])
set(gca,'XTick',[0:4:24]);
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
xlim([0 24])
ylabel('[kW]');
set(gca,'XTick',[0:4:24]);
xlabel('Time [h]')
box on


%%
figure
%plot((t_start:t_start+Tsim-1)/60,100*datalog.SOC(2,1:end), ':','linewidth',2.5)
hold on
plot((t_start:t_start+Tsim-1)/60,datalog.OPFsim.p_G_sim(batt_nodes(2),1:end),'-r','linewidth',2.5)
hold on
plot((t_start:t_start+Tsim-1)/60,datalog.MPC.P_ref(batt_nodes(2),1:end), '--b','linewidth',2.5)
%title('Batt 2')
grid on 
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
%legend({'Real power','Reference power'},'Orientation','horizontal','location','north');
xlim([0 24])
set(gca,'XTick',[0:4:24]);
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
xlim([0 24])
ylabel('[kW]');
set(gca,'XTick',[0:4:24]);
xlabel('Time [h]')
box on

%%
figure
%plot((t_start:t_start+Tsim-1)/60,datalog.OPFsim.p_G_sim(batt_nodes(3),1:end),'-r','linewidth',2.5)
hold on
plot((t_start:t_start+Tsim-1)/60,datalog.OPFsim.p_G_sim(batt_nodes(3),1:end),'-r','linewidth',2.5)
hold on
plot((t_start:t_start+Tsim-1)/60,datalog.MPC.P_ref(batt_nodes(3),1:end), '--b','linewidth',2.5)
%title('Batt 3')
grid on 
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
%legend({'Real power','Reference power'},'Orientation','horizontal','location','north');
xlim([0 24])
set(gca,'XTick',[0:4:24]);
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
xlim([0 24])
ylabel('[kW]');
xlabel('Time [h]')
box on

%%

figure
plot((t_start:t_start+Tsim-1)/60,datalog.SOC(1,1:end), ':','linewidth',2.5)
hold on
plot((t_start:t_start+Tsim-1)/60,datalog.SOC(2,1:end), '--','linewidth',2.5)
hold on
plot((t_start:t_start+Tsim-1)/60,datalog.SOC(3,1:end), '-','linewidth',2.5)
grid on 
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
%legend({'Batt 1','Batt 2','Batt 3'},'Orientation','horizontal','location','south');
xlim([0 24])
set(gca,'XTick',[0:4:24]);
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
xlim([0 24])
ylim([0 1])
%ylabel('[kW]');
xlabel('Time [h]')
box on

%%
figure
plot((t_start:t_start+Tsim-1)/60,datalog.OPFsim.p_G_sim(gen_nodes(1),1:end),'-r','linewidth',2.5)
hold on
plot((t_start:t_start+Tsim-1)/60,datalog.MPC.P_ref(gen_nodes(1),1:end),'--b','linewidth',2.5)
hold on
%plot((t_start:t_start+Tsim-1)/60,datalog.OPF.Pg(gen_nodes(:),1:end),'-r')
hold on
%title('Dispatchable generator 1')
grid on 
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
xlim([0 24])
ylabel('[kW]');
ylim([0 60])
set(gca,'XTick',[0:4:24]);
xlabel('Time [h]')
box on

figure
plot((t_start:t_start+Tsim-1)/60,datalog.OPFsim.p_G_sim(gen_nodes(2),1:end),'-r','linewidth',2.5)
hold on
plot((t_start:t_start+Tsim-1)/60,datalog.MPC.P_ref(gen_nodes(2),1:end),'--b','linewidth',2.5)
hold on
%plot((t_start:t_start+Tsim-1)/60,datalog.OPF.Pg(gen_nodes(:),1:end),'-r')
hold on
%title('Dispatchable generator 2')
grid on 
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
xlim([0 24])
ylim([0 60])
ylabel('[kW]');
set(gca,'XTick',[0:4:24]); 
xlabel('Time [h]')
box on


%%


figure
plot((t_start:t_start+Tsim-1)/60, 90*ones(size(datalog.OPFsim.v_sim(2:end,:),1),size(datalog.OPFsim.v_sim(2:end,:),2)),'--k','linewidth',1)
hold on
plot((t_start:t_start+Tsim-1)/60, datalog.OPFsim.v_sim(2:end,:),'linewidth',1)
hold on
plot((t_start:t_start+Tsim-1)/60, 110*ones(size(datalog.OPFsim.v_sim(2:end,:),1),size(datalog.OPFsim.v_sim(2:end,:),2)),'--k','linewidth',1)
grid on 
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
xlim([0 24])
ylim([80 120])
ylabel('[V]')
set(gca,'XTick',[0:4:24]);
xlabel('Time [h]')
box on

%%

load data_sim_LoadMeasure
 

sim_init = 554;
sim_centr = 555;
sim_end  = 556;

sourcenodes_mod = 2:6;
figure
plot(  t_plot(:,(sim_init)*(length(tr))+1  :(sim_end)*(length(tr)))/3600 ,   y_plot(sourcenodes_mod,(sim_init-1)*(length(tr))+1 :(sim_end-1)*(length(tr)) ),'linewidth',2.5)
hold on
plot(  t_plot(:,(sim_init)*(length(tr))+1  :(sim_end)*(length(tr)))/3600 ,  kron(datalog.OPF.v_opf(sourcenodes_mod,sim_init:sim_end-1),ones(1,length(tr))),':k','linewidth',2.5)
xlim( [sim_centr/60-2/(60*60), sim_centr/60+2/(60*60)])
ylim([min(min(datalog.OPF.v_opf(sourcenodes_mod,sim_init:sim_end-1)))-0.5, max(max(datalog.OPF.v_opf(sourcenodes_mod,sim_init:sim_end-1))+0.5)])
grid on
ylabel('[V]')
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
xlabel('Time [h]')
box on

figure
plot(  t_plot(:,(sim_init)*(length(tr))+1  :(sim_end)*(length(tr)))/3600 ,   y_plot(loadnodes,(sim_init-1)*(length(tr))+1 :(sim_end-1)*(length(tr)) ),'linewidth',2.5)
hold on
%plot(  t_plot(:,(569)*(length(tr))+1  :(571)*(length(tr)))/3600 ,  kron(datalog.OPF.v_opf(2,569:570),ones(1,length(tr))),'--k','linewidth',2.5)
xlim( [sim_centr/60-2/(60*60), sim_centr/60+2/(60*60)])
%ylim([min(datalog.OPF.v_opf(2,sim_init:sim_end-1))-0.5, max(datalog.OPF.v_opf(2,sim_init:sim_end-1))+0.5])
grid on
ylabel('[V]')
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
xlabel('Time [h]')
box on

