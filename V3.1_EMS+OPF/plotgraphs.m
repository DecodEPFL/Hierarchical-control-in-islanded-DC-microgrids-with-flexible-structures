close all;

%%
figure(1);
plot(t_plot,y_plot(:,1:N),'linewidth',2.5);
set(gca, 'FontName', 'Times New Roman')
h=legend('$V_1$','$V_2$','$V_3$','$V_4$','$V_5$','$V_6$','Location','north');
set(h, 'Interpreter', 'latex')
xlabel('Time (sec)')
ylabel('PCC voltage (V)')
set(gca,'fontsize',15)
set(gca,'XTick',[0:2*t_sim:Tot_sim]);
axis([0,Tot_sim,min(min(y_plot(:,1:N)))-0.5,max(max(y_plot(:,1:N)))+1])
grid on
hold off plots
%%
% figure(2);
% plot(t_plot,y_plot(:,N+1:N+n),'linewidth',1.5);
% set(gca, 'FontName', 'Times New Roman')
% legend('I_{t1}','I_{t2}','Location','southeast')
% xlabel('Time (seconds)')
% ylabel('(A)')
% grid on
% hold off plots

%%
P_gen = y_plot(:,N+1:N+n).*y_plot(:,1:n) + y_plot(:,N+1:N+n).* y_plot(:,N+1:N+n)*R_t;

figure(3);
h1=plot(t_plot,P_gen,'linewidth',2);
set(h1, {'color'}, {[0.7500         0    0.7500]; [0 0.5 0]});
hold on
h2=plot(t_plot,P_G_tot_ref_plot(1,:),'--r','linewidth',2.5);
hold on
h3=plot(t_plot,P_G_tot_ref_plot(2,:),'--r','linewidth',2.5, 'HandleVisibility','off');
hold on
h4=plot(t_plot,kron(p_G_ref(1,1)',ones(length(t_plot),1)),':b','linewidth',1.5);
hold on
h5=plot(t_plot,kron(p_G_ref(2,1)',ones(length(t_plot),1)),':b','linewidth',1.5, 'HandleVisibility','off');
set(gca, 'FontName', 'Times New Roman')
h=legend('$P_{t1}$','$P_{t2}$','${P}_{G}^*$','$\bar{P}_{G}$','Location','south','Orientation','horizontal');
set(h, 'Interpreter', 'latex');
xlabel('Time (seconds)')
ylabel('(W)')
set(gca,'fontsize',15)
set(gca,'XTick',[0:2*t_sim:Tot_sim]);
axis([0,Tot_sim,min(min(P_gen))-20,max(max(P_gen))+20])
grid on
hold off plots
%%
% figure(4);
% plot(t_plot,y_plot(:,N+2*n+1:N+2*n+M),'linewidth',1.5);
% set(gca, 'FontName', 'Times New Roman')
% legend('I_{1}','I_{2}','I_{3}','I_{4}','I_{5}','I_{6}','I_{7}','Location','southeast')
% xlabel('Time (seconds)')
% ylabel('(A)')
% grid on
% hold off plots

%%
figure(5);
subplot 211
plot(t_plot,P_L_const_trend(loadnodes,:),'linewidth',3);
set(gca, 'FontName', 'Times New Roman')
ylabel('(W)')
axis([ 0 t_plot(end) 0 80])
h=legend('$P_{L3}$','$P_{L4}$','$P_{L5}$','$P_{L6}$','Location','north','Orientation','horizontal');
set(h, 'Interpreter', 'latex');
xlabel('Time (seconds)')
grid on
set(gca,'fontsize',15)
set(gca,'XTick',[0,t_sim:2*t_sim:Tot_sim]);
axis([0,Tot_sim,min(min(P_L_const_trend))-20,max(max(P_L_const_trend))+30])

grid on
subplot 212
plot(t_plot,I_L_const_trend(loadnodes,:),'linewidth',3);
set(gca, 'FontName', 'Times New Roman')
ylabel('(A)')
axis([ 0 t_plot(end) 0 10])
% legend('Pload-1','Pload-2','Pload-3','Pload-4','Pload-5','Pload-6','Location','southeast')
xlabel('Time (seconds)')
grid on
set(gca,'fontsize',15)
set(gca,'XTick',[0,t_sim:2*t_sim:Tot_sim]);
axis([0,Tot_sim,min(min(I_L_const_trend))-1,max(max(I_L_const_trend))+4])
h=legend('$\bar{I}_{L3}$','$\bar{I}_{L4}$','$\bar{I}_{L5}$','$\bar{I}_{L6}$','Location','north','Orientation','horizontal');
set(h, 'Interpreter', 'latex');
hold off plots