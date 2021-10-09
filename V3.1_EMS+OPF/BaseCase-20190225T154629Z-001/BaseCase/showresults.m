%% Results

Initialization

load datalog

prompt1=('Initial instant to plot?           ');

a=input(prompt1);

prompt2=('Final instant to plot?             ');

Tt=input(prompt2);

% Parameters loading and definition

TestFacility_RealData

n_nodes    =     netdata.n_nodes;
nb         =     length(gen_batt);                                         % Number of batteries definiton
t          =     size(datalog.X,2);                                        % Number of simulated steps

% Fundamental variables to plot are extracted from the datalog

w          =     datalog.X(end-1-nb,1:t);

Vp=zeros(n_nodes,t);
Vp(1:n_nodes,:) = datalog.X(1:n_nodes,1:t);

angleV = datalog.X(2:n_nodes,1:t);

v = datalog.X(end-nb,1:t);

Soc = datalog.X(end-nb+1:end,1:t);


% The figures are now reported

figure
plot((a/120:1/120:Tt/120),(Soc(1,a:Tt)),'b',(a/120:1/120:Tt/120),(Soc(2,a:Tt)),'g',(a/120:1/120:Tt/120),(Soc(3,a:Tt)),'r','LineWidth',1)
legend('Battery 1 (Node 5)','Battery 2 (Node 9)','Battery 3 (Node 10)','Location','NorthEast')
hold on
plot((a/120:1/120:Tt/120),20*ones(1,Tt-a+1),'--k',(a/120:1/120:Tt/120),80*ones(1,Tt-a+1),'--k','LineWidth',0.5)
title('SOCs')
ylabel('States of Charge [%] ')
xlabel(' Time [hours] ')
ylim([0, 100])
xlim([a/120, Tt/120])
grid

figure
plot(((a/120:1/120:Tt/120)),w(a:Tt)/2/pi,'LineWidth',1)
hold on
plot((a/120:1/120:Tt/120),49*ones(1,Tt-a+1),'--k',(a/120:1/120:Tt/120),51*ones(1,Tt-a+1),'--k','LineWidth',0.5)
title('Frequency')
ylabel(' Network Frequency [Hz] ')
xlabel(' Time [hours] ')
ylim([48, 52])
xlim([a/120, Tt/120])
grid

figure
plot(((a/120:1/120:Tt/120)),Vp(1,(a:Tt)),((a/120:1/120:Tt/120)),Vp(2,(a:Tt)),((a/120:1/120:Tt/120)),Vp(3,(a:Tt)),((a/120:1/120:Tt/120)),Vp(4,(a:Tt)),((a/120:1/120:Tt/120)),Vp(5,a:Tt),((a/120:1/120:Tt/120)),Vp(6,a:Tt),((a/120:1/120:Tt/120)),Vp(7,a:Tt),((a/120:1/120:Tt/120)),Vp(8,a:Tt),((a/120:1/120:Tt/120)),Vp(9,a:Tt),((a/120:1/120:Tt/120)),Vp(10,a:Tt),((a/120:1/120:Tt/120)),Vp(11,a:Tt),((a/120:1/120:Tt/120)),Vp(12,a:Tt),((a/120:1/120:Tt/120)),Vp(13,a:Tt),'LineWidth',1);
hold on
plot((a/120:1/120:Tt/120),360*ones(1,Tt-a+1),'--k',(a/120:1/120:Tt/120),440*ones(1,Tt-a+1),'--k','LineWidth',0.5)
ylim([350 , 450]);
title('Voltages')
legend('Nodal Voltage 1','Nodal Voltage 2','Nodal Voltage 3','Nodal Voltage 4','Nodal Voltage 5','Nodal Voltage 6','Nodal Voltage 7','Nodal Voltage 8','Nodal Voltage 9','Nodal Voltage 10','Nodal Voltage 11','Nodal Voltage 12','Nodal Voltage 13','Location','NorthWest')
ylabel('Nodal Voltages [V] ')
xlabel(' Time [hours] ')
xlim([a/120, Tt/120])
grid


figure
plot(((a/120:1/120:Tt/120)),datalog.U(1,a:Tt)/1000,'b',((a/120:1/120:Tt/120)),datalog.U(2,a:Tt)/1000,'g',((a/120:1/120:Tt/120)),datalog.U(3,a:Tt)/1000,'r',((a/120:1/120:Tt/120)),datalog.U(4,a:Tt)/1000,'m','LineWidth',1)
hold on
plot(((a/120:1/120:Tt/120)),real(Pgen(a:Tt,5))/1000,'Color',[1 0.84 0],'LineWidth',1)
hold on
plot(((a/120:1/120:Tt/120)),real(Pgen(a:Tt,6))/1000,'Color',[0.3 0.75 0.93],'LineWidth',1)
title(' Reference Powers: Active Power ' )
legend('Battery 1 (Node 5)','Battery 2 (Node 9)','Battery 3 (Node 10)','Rot. Generator (Node 11)','Solar Panel (Node 12)','Wind Turbine (Node 13)','Location','NorthEast')
ylabel(' Active Power [kW] ')
xlabel(' Time [hours] ')
xlim([a/120, Tt/120])
grid

figure
plot(((a/120:1/120:Tt/120)),datalog.U(5,a:Tt)/1000,'b',((a/120:1/120:Tt/120)),datalog.U(6,a:Tt)/1000,'g',((a/120:1/120:Tt/120)),datalog.U(7,a:Tt)/1000,'r',((a/120:1/120:Tt/120)),datalog.U(8,a:Tt)/1000,'m','LineWidth',1)
hold on
plot(((a/120:1/120:Tt/120)),imag(Pgen(a:Tt,5))/1000,'Color',[1 0.84 0],'LineWidth',1)
hold on
plot(((a/120:1/120:Tt/120)),imag(Pgen(a:Tt,6))/1000,'Color',[0.3 0.75 0.93],'LineWidth',1)
legend('Battery 1 (Node 5)','Battery 2 (Node 9)','Battery 3 (Node 10)','Rot. Generator (Node 11)','Solar Panel (Node 12)','Wind Turbine (Node 13)','Location','NorthEast')
title(' Reference Powers: Reactive Power ' )
ylabel(' Reactive Power [kVar] ')
xlabel(' Time [hours] ')
xlim([a/120, Tt/120])
grid

figure
plot(((a/120:1/120:Tt/120)),real(datalog.Sgen(1,a:Tt))/1000,'b',((a/120:1/120:Tt/120)),real(datalog.Sgen(2,a:Tt))/1000,'g',((a/120:1/120:Tt/120)),real(datalog.Sgen(3,a:Tt))/1000,'r',((a/120:1/120:Tt/120)),real(datalog.Sgen(4,a:Tt))/1000,'m','LineWidth',1)
hold on
plot(((a/120:1/120:Tt/120)),real(datalog.Sgen(5,a:Tt))/1000,'Color',[1 0.84 0],'LineWidth',1)
hold on
plot((a/120:1/120:Tt/120),real(datalog.Sgen(6,a:Tt))/1000,'Color',[0.3 0.75 0.93],'LineWidth',1)
legend('Battery 1 (Node 5)','Battery 2 (Node 9)','Battery 3 (Node 10)','Rot. Generator (Node 11)','Solar Panel (Node 12)','Wind Turbine (Node 13)','Location','NorthEast')
title(' Generated Powers: Active Power ' )
ylabel(' Active Power [kW] ')
xlabel(' Time [hours] ')
xlim([a/120, Tt/120])
grid


figure
plot(((a/120:1/120:Tt/120)),imag(datalog.Sgen(1,a:Tt))/1000,'b',((a/120:1/120:Tt/120)),imag(datalog.Sgen(2,a:Tt))/1000,'g',((a/120:1/120:Tt/120)),imag(datalog.Sgen(3,a:Tt))/1000,'r',(a/120:1/120:Tt/120),imag(datalog.Sgen(4,a:Tt))/1000,'m','LineWidth',1)
hold on
plot((a/120:1/120:Tt/120),imag(datalog.Sgen(5,a:Tt))/1000,'Color',[1 0.84 0],'LineWidth',1)
hold on
plot((a/120:1/120:Tt/120),imag(datalog.Sgen(6,a:Tt))/1000,'Color',[0.3 0.75 0.93],'LineWidth',1)
legend('Battery 1 (Node 5)','Battery 2 (Node 9)','Battery 3 (Node 10)','Rot. Generator (Node 11)','Solar Panel (Node 12)','Wind Turbine (Node 13)','Location','NorthEast')
title(' Generated Powers: Reactive Power ' )
ylabel(' Reactive Power [kVar] ')
xlabel(' Time [hours] ')
xlim([a/120, Tt/120])
grid

figure
plot((a/120:1/120:Tt/120),real(Pload(a:Tt,1))/1000,'b',(a/120:1/120:Tt/120),real(Pload(a:Tt,2))/1000,'g',(a/120:1/120:Tt/120),real(Pload(a:Tt,3))/1000,'r',(a/120:1/120:Tt/120),real(Pload(a:Tt,4))/1000,'y',(a/120:1/120:Tt/120),real(Pload(a:Tt,5))/1000,'m')
title('Load Profile: Active Power');
ylim([-6 , +50]);
ylabel(' Active Power [kW] ')
xlabel(' Time [hours] ')
xlim([a/120, Tt/120])
grid

figure
plot((a/120:1/120:Tt/120),imag(Pload(a:Tt,1))/1000,'b',(a/120:1/120:Tt/120),imag(Pload(a:Tt,2))/1000,'g',(a/120:1/120:Tt/120),imag(Pload(a:Tt,3))/1000,'r',(a/120:1/120:Tt/120),imag(Pload(a:Tt,4))/1000,'y',(a/120:1/120:Tt/120),imag(Pload(a:Tt,5))/1000,'m')
title('Load Profile: Reactive Power');
ylim([-6 , +50]);
ylabel(' Reactive Power [kVar] ')
xlabel(' Time [hours] ')
xlim([a/120, Tt/120])
grid

figure
plot((a/120:1/120:Tt/120),real(datalog.Sload(1,a:Tt))/1000,'b',(a/120:1/120:Tt/120),real(datalog.Sload(2,a:Tt))/1000,'g',(a/120:1/120:Tt/120),real(datalog.Sload(3,a:Tt))/1000,'r',(a/120:1/120:Tt/120),real(datalog.Sload(4,a:Tt))/1000,'y',(a/120:1/120:Tt/120),real(datalog.Sload(5,a:Tt))/1000,'m')
title('Load Power: Active Power');
ylim([-6 , +50]);
ylabel(' Active Power [kW] ')
xlabel(' Time [hours] ')
xlim([a/120, Tt/120])
grid

figure
plot((a/120:1/120:Tt/120),imag(datalog.Sload(1,a:Tt))/1000,'b',(a/120:1/120:Tt/120),imag(datalog.Sload(2,a:Tt))/1000,'g',(a/120:1/120:Tt/120),imag(datalog.Sload(3,a:Tt))/1000,'r',(a/120:1/120:Tt/120),imag(datalog.Sload(4,a:Tt))/1000,'y',(a/120:1/120:Tt/120),imag(datalog.Sload(5,a:Tt))/1000,'m')
title('Load Power: Reactive Power');
ylim([-6 , +50]);
ylabel(' Reactive Power [kVar] ')
xlabel(' Time [hours] ')
xlim([a/120, Tt/120])
grid



