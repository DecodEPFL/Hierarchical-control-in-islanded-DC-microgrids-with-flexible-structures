%% LOAD FLOW
Dati_rete

generatori_e_carichi
results = connected_PF(netdata,RegG,RegL);

Vlf = abs(results.V);
Vlf_V=Vlf*Vrif;
angleVlf = angle(results.V);
angleVlf_deg = angleVlf*180/pi;

Plf = real(results.S);
Qlf = imag(results.S);
Plf_kW = real(results.S)*Srif/1000;
Qlf_kVar = imag(results.S)*Srif/1000;

%% Dati i risultati del load flow, testiamo lo Jacobiano

[J,Pt,Qt,Y] = RealJacobian(Vlf,angleVlf,netdata)

%% Test Giacomino
RegL(1).P=-1250e3/Srif;
RegL(2).P=-500e3/Srif;
RegG(1).P = 500e3/Srif;
RegG(2).P = 500e3/Srif;

results2=connected_PF(netdata,RegG,RegL);

Vlf2 = abs(results2.V);
Vlf2_V = Vlf2*Vrif;
angleVlf2 = angle(results2.V);
angleVlf2_deg = angleVlf2*180/pi;

Plf2 = real(results2.S);
Qlf2 = imag(results2.S);
Plf2_kW = real(results2.S)*Srif/1000;
Qlf2_kVar = imag(results2.S)*Srif/1000;

DeltaP = (Plf2(2:end)-Plf(2:end));
DeltaQ = (Qlf2(2:end)-Qlf(2:end));

Jinv = inv(J);

new_deltax=Jinv*[DeltaP;DeltaQ];

new_deltaV=new_deltax(1:5);
new_deltaangleV=new_deltax(6:end);

V2 = +[0;new_deltaV]+Vlf;
V2_V = V2*Vrif;

angle2 = [0;new_deltaangleV]+angleVlf;
angle2_deg = angle2*180/pi;

dV_J=V2_V-Vlf_V
dV_lf=Vlf2_V-Vlf_V

e_V = (Vlf2_V-V2_V)
e_V_percent = e_V./dV_lf*100
e_delta = (angleVlf2 - angle2)
e_delta_deg = (angleVlf2_deg - angle2_deg)

%% Potenza nodo di saldo

dPs_dV2 = abs(Y(1,2))*cos(-angleVlf(2)-angle(Y(1,2)));

dPs_dangle2 = Vlf(2)*abs(Y(1,2))*sin(-angleVlf(2)-angle(Y(1,2)));

Ps = Plf(1) + dPs_dV2*new_deltaV(1) + dPs_dangle2*new_deltaangleV(1);

Ps_new = Ps*Srif;