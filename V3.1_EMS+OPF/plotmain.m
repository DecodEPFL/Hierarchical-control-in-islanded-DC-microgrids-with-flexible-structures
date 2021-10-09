%%
close all

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
plot((t_start:t_start+Tsim-1)/60,realProfiles.L(1).I, 'linewidth',1.5)
hold on
plot((t_start:t_start+Tsim-1)/60,realProfiles.L(4).I, 'linewidth',1.5)
hold on
plot((t_start:t_start+Tsim-1)/60,realProfiles.L(7).I, 'linewidth',1.5)
grid on 
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
xlim([0 24])
set(gca,'XTick',[0:4:24]);

figure
plot((t_start:t_start+Tsim-1)/60,realProfiles.L(1).P, 'linewidth',1.5)
hold on
plot((t_start:t_start+Tsim-1)/60,realProfiles.L(4).P, 'linewidth',1.5)
hold on
plot((t_start:t_start+Tsim-1)/60,realProfiles.L(7).P, 'linewidth',1.5)
grid on 
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
xlim([0 24])
set(gca,'XTick',[0:4:24]);
%%

figure
plot((t_start:t_start+Tsim-1)/60, datalog.OPF.v_opf(:,:))
hold on
plot((t_start:t_start+Tsim-1)/60, datalog.OPFsim.v_sim(:,:))
grid on 
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
xlim([0 24])
set(gca,'XTick',[0:4:24]);
%%

figure
%plot((t_start:t_start+Tsim-1)/60,P_nom(2,t_start : (t_start+Tsim)-1), 'linewidth',1.5)
%hold on
plot((t_start:t_start+Tsim-1)/60,P_real(7, t_start : (t_start+Tsim)-1), 'linewidth',1.5)
hold on
plot((t_start:15:t_start+Tsim-1)/60,PL_nom_s(7,1 : 96), '--k', 'linewidth',2.5)
title('PL nom2')
grid on 
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
xlim([0 24])
set(gca,'XTick',[0:4:24]);


%%
figure
plot((t_start:t_start+Tsim-1)/60,-P_nom(1,t_start : (t_start+Tsim)-1), 'linewidth',1.5)
hold on
plot((t_start:t_start+Tsim-1)/60, datalog.MPC.P_ref(end, t_start : (t_start+Tsim)-1)  , '--k', 'linewidth',1.5)
hold on
plot((t_start*MPCpar.t_s/60: MPCpar.t_s/60 :(t_start+Tsim-1)/60 ), -PL_nom_s(1, (t_start-1)*MPCpar.t_s + 1 : Tsim/MPCpar.t_s ), ':k', 'linewidth',1.5)
hold on
plot((t_start:t_start+Tsim-1)/60,datalog.OPF.Pg(pvnodes,1:end),'-r')
% legend('Nominal power','Actual power');
title('PV system')
grid on 
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
xlim([0 24])
set(gca,'XTick',[0:4:24]);

%%
figure
plot((t_start:t_start+Tsim-1)/60,100*datalog.SOC(1,1:end), 'linewidth',1.5)
hold on
plot((t_start:t_start+Tsim-1)/60,datalog.MPC.P_ref(batt_nodes(1),1:end),'-b', 'linewidth',1)
hold on
plot((t_start:t_start+Tsim-1)/60,datalog.OPFsim.p_G_sim(batt_nodes(1),1:end),'-r','linewidth',1)
title('Batt 1')
grid on 
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
legend('SOC1','Pb1');
% xlim([0 24])
% set(gca,'XTick',[0:4:24]);

%%
figure
plot((t_start:t_start+Tsim-1)/60,100*datalog.SOC(2,1:end), 'linewidth',1.5)
hold on
plot((t_start:t_start+Tsim-1)/60,datalog.MPC.P_ref(batt_nodes(2),1:end), '-b')
hold on
plot((t_start:t_start+Tsim-1)/60,datalog.OPFsim.p_G_sim(batt_nodes(2),1:end),'-r')
title('Batt 2')
grid on 
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
legend('SOC2','Pb2');
% xlim([0 24])
% set(gca,'XTick',[0:4:24]);

%%
figure
plot((t_start:t_start+Tsim-1)/60,100*datalog.SOC(3,1:end), 'linewidth',1.5)
hold on
plot((t_start:t_start+Tsim-1)/60,datalog.MPC.P_ref(batt_nodes(3),1:end), '-b')
hold on
plot((t_start:t_start+Tsim-1)/60,datalog.OPFsim.p_G_sim(batt_nodes(3),1:end),'-r')
title('Batt 3')
grid on 
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
legend('SOC3','Pb3');
% xlim([0 24])
% set(gca,'XTick',[0:4:24]);

%%
figure
plot((t_start:t_start+Tsim-1)/60,datalog.MPC.P_ref(gen_nodes,1:end), 'linewidth',1)
hold on
plot((t_start:t_start+Tsim-1)/60,datalog.OPFsim.p_G_sim(gen_nodes,1:end),'-r')
title('P - gen')
grid on 
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
% xlim([0 24])
% set(gca,'XTick',[0:4:24]);
% 
%%
figure
plot((t_start:t_start+Tsim-1)/60,-P_nom(1,t_start : (t_start+Tsim)-1), 'linewidth',1.5)
hold on
plot((t_start:t_start+Tsim-1)/60, datalog.MPC.P_ref(end, t_start : (t_start+Tsim)-1)  , '--k', 'linewidth',1.5)
hold on
plot((t_start*MPCpar.t_s/60: MPCpar.t_s/60 :(t_start+Tsim-1)/60 ), -PL_nom_s(1, (t_start-1)*MPCpar.t_s + 1 : Tsim/MPCpar.t_s ), ':k', 'linewidth',1.5)
hold on
plot((t_start:t_start+Tsim-1)/60,datalog.OPFsim.p_G_sim(pvnodes,1:end),'-r')
% legend('Nominal power','Actual power');
title('PV system')
grid on 
set(gca,'fontsize',15, 'FontName', 'Times New Roman')
xlim([0 24])
set(gca,'XTick',[0:4:24]);

 %%
% 
% figure
% plot((t_start:t_start+Tsim-1)/60,P_nom(3,t_start : (t_start+Tsim)-1), 'linewidth',1.5)
% title('PL nom3')
% grid on 
% set(gca,'fontsize',15, 'FontName', 'Times New Roman')
% xlim([0 24])
% set(gca,'XTick',[0:4:24]);
% 
% figure
% plot((t_start:t_start+Tsim-1)/60,P_nom(4,t_start : (t_start+Tsim)-1), 'linewidth',1.5)
% title('PL nom4')
% grid on 
% set(gca,'fontsize',15, 'FontName', 'Times New Roman')
% xlim([0 24])
% set(gca,'XTick',[0:4:24]);