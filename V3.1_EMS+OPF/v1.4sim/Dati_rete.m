
% dati di rete, generatori e carichi
netdata = struct;
%riferimenti pu

Vrif=20e3;
Srif=1e6;
Zrif=Vrif^2/Srif;

netdata.Vrif = Vrif;
netdata.Srif = Srif;
netdata.Zrif = Zrif;

% ni - vettore nodi inizio linea
% nf - vettore nodi fine linea

netdata.n_nodes  =    6;

netdata.ni       =    [1,  2,   3,  4,  4];

netdata.nf       =    [2,  3,   4,  5,  6];

%% dati rete

f=50;

rho_XoverR = 2.5;

R12=5;
L12=R12*rho_XoverR/(2*pi*f);

R23=5;
L23=R23*rho_XoverR/(2*pi*f);

R34=5;
L34=R34*rho_XoverR/(2*pi*f);

R45=5;
L45=R45*rho_XoverR/(2*pi*f);

R46=5;
L46=R46*rho_XoverR/(2*pi*f);

netdata.nodi_ctrl_gen        =      [5, 6];
netdata.nodi_batt_pq         =      [5, 6];
netdata.nodi_load            =      [2, 3];
netdata.nodi_non_crtl_gen    =           5;
netdata.nodi_gen             =   [5, 5, 6];    

% nodo di saldo e sua tensione [V]
netdata.ns=1;
netdata.Vs=20e3/Vrif;
netdata.Vnom =20e3/Vrif;

% frequenza nominale [rad/s]
netdata.wn=2*pi*50;

% Parametri circuitali delle linee
% R  - resistenza  [Ohm]
% L  - induttanza  [H]
% C  - capacità    [F]
% G  - conduttanza [Siemens]
% T  - elastanza   [1/H]
% k  - rapporto di trasformazione [V/V]
% 
% netdata.r            =    [R12;        R23;            R34;         R45;            R46]/Zrif;   % pu/km
%  
% netdata.l            =    [L12;        L23;            L34;         L45;            L46]/Zrif; % pu/km
%  
% netdata.linelength   =    [2;            1;              3;           1;              1];        % km


netdata.r            =    [0.64;        1;            0.5;         1;            1.5]/Zrif;   % pu/km

netdata.l            =    [0.0815;  2/(2*pi*50);  1/(2*pi*50);  1.5/(2*pi*50);   2/(2*pi*50)]/Zrif; % pu/km

netdata.linelength   =    [1;           1;              1;           1;              1];        % km



netdata.R            =    (netdata.r).*netdata.linelength;       % pu

netdata.L            =    (netdata.l).*netdata.linelength;       % pu

netdata.C            =    zeros(size(netdata.R));

netdata.G            =    zeros(size(netdata.R));
 
netdata.T            =    zeros(size(netdata.R));

netdata.k            =    ones(size(netdata.R));                                    % They are all set to 1 since there are no trasformers in the considered microgrid


% Dati per calcolo perdite per modello linee a pi-greco (Vedi libro Tironi)

alpha1_12      =  abs(1 );% +   (netdata.R(1) + 1i*netdata.L(1)*(2*pi*50)) * (netdata.G(1) + 1/(1i*netdata.C(1)*(2*pi*50)))/2   );
zeta1_12       =  abs((netdata.R(1) + 1i*netdata.L(1)*(2*pi*50)));
phiz1_12        =  angle((netdata.R(1) + 1i*netdata.L(1)*(2*pi*50)));
phialpha1_12   =  angle(1 );% +  (netdata.R(1) + 1i*netdata.L(1)*(2*pi*50)) * (netdata.G(1) + 1/(1i*netdata.C(1)*(2*pi*50)))/2   );

alpha1_23      =  abs(1 );%   (netdata.R(2) + 1i*netdata.L(2)*(2*pi*50)) * (netdata.G(2) + 1/(1i*netdata.C(2)*(2*pi*50)))/2   );
zeta1_23       =  abs((netdata.R(2) + 1i*netdata.L(2)*(2*pi*50)));
phiz1_23        =  angle((netdata.R(2) + 1i*netdata.L(2)*(2*pi*50)));
phialpha1_23   =  angle(1 );%+   (netdata.R(2) + 1i*netdata.L(2)*(2*pi*50)) * (netdata.G(2) + 1/(1i*netdata.C(2)*(2*pi*50)))/2   );

alpha1_34      =  abs(1 );%+   (netdata.R(3) + 1i*netdata.L(3)*(2*pi*50)) * (netdata.G(3) + 1/(1i*netdata.C(3)*(2*pi*50)))/2   );
zeta1_34       =  abs((netdata.R(3) + 1i*netdata.L(3)*(2*pi*50)));
phiz1_34        =  angle((netdata.R(3) + 1i*netdata.L(3)*(2*pi*50)));
phialpha1_34   =  angle(1 );%+   (netdata.R(3) + 1i*netdata.L(3)*(2*pi*50)) * (netdata.G(3) + 1/(1i*netdata.C(3)*(2*pi*50)))/2   );

alpha1_45      =  abs(1 );%+   (netdata.R(4) + 1i*netdata.L(4)*(2*pi*50)) * (netdata.G(4) + 1/(1i*netdata.C(4)*(2*pi*50)))/2   );
zeta1_45       =  abs((netdata.R(4) + 1i*netdata.L(4)*(2*pi*50)));
phiz1_45        =  angle((netdata.R(4) + 1i*netdata.L(4)*(2*pi*50)));
phialpha1_45   =  angle(1 );%+   (netdata.R(4) + 1i*netdata.L(4)*(2*pi*50)) * (netdata.G(4) + 1/(1i*netdata.C(4)*(2*pi*50)))/2   );

alpha1_46      =  abs(1);% +   (netdata.R(5) + 1i*netdata.L(5)*(2*pi*50)) * (netdata.G(5) + 1/(1i*netdata.C(5)*(2*pi*50)))/2   );
zeta1_46       =  abs((netdata.R(5) + 1i*netdata.L(5)*(2*pi*50)));
phiz1_46        =  angle((netdata.R(5) + 1i*netdata.L(5)*(2*pi*50)));
phialpha1_46   =  angle(1 );%+   (netdata.R(5) + 1i*netdata.L(5)*(2*pi*50)) * (netdata.G(5) + 1/(1i*netdata.C(5)*(2*pi*50)))/2   );


% Caratteristiche rete

netdata.b12       =   (alpha1_12/zeta1_12)*cos(phiz1_12 - phialpha1_12);
netdata.gamma12   =   (2/zeta1_12)*cos(phiz1_12);

netdata.b23       =   (alpha1_23/zeta1_23)*cos(phiz1_23 - phialpha1_23);
netdata.gamma23   =   (2/zeta1_23)*cos(phiz1_23);

netdata.b34       =   (alpha1_34/zeta1_34)*cos(phiz1_34 - phialpha1_34);
netdata.gamma34   =   (2/zeta1_34)*cos(phiz1_34);

netdata.b45       =   (alpha1_45/zeta1_45)*cos(phiz1_45 - phialpha1_45);
netdata.gamma45   =   (2/zeta1_45)*cos(phiz1_45);

netdata.b46       =   (alpha1_46/zeta1_46)*cos(phiz1_46 - phialpha1_46);
netdata.gamma46   =   (2/zeta1_46)*cos(phiz1_46);

