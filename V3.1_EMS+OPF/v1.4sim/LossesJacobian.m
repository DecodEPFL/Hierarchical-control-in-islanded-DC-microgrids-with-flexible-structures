function  [Jloss, Ploss] = LossesJacobian(V, angleV, netdata)


n = netdata.n_nodes;


Jloss = zeros(n-1,2*n-2);

%Ploss12
Jloss(1,1)   =       2*netdata.b12*V(2) - netdata.gamma12*V(1)*cos(angleV(1) - angleV(2));
Jloss(1,n)   =     - netdata.gamma12*V(1)*V(2)*sin(angleV(1) - angleV(2));

%Ploss23
Jloss(2,1)     =     2*netdata.b23*V(2) - netdata.gamma23*V(3)*cos(angleV(2) - angleV(3));
Jloss(2,2)     =     2*netdata.b23*V(3) - netdata.gamma23*V(2)*cos(angleV(2) - angleV(3));
Jloss(2,n)     =     netdata.gamma23*V(2)*V(3)*sin(angleV(2) - angleV(3));
Jloss(2,n+1)   =    - Jloss(2,n);

%Ploss34
Jloss(3,2)     =     2*netdata.b34*V(3) - netdata.gamma34*V(4)*cos(angleV(3) - angleV(4));
Jloss(3,3)     =     2*netdata.b34*V(4) - netdata.gamma34*V(3)*cos(angleV(3) - angleV(4));
Jloss(3,n+1)   =     netdata.gamma34*V(4)*V(3)*sin(angleV(3) - angleV(4));
Jloss(3,n+2)   =    -Jloss(3,n+1);

%Ploss45
Jloss(4,3)     =     2*netdata.b45*V(4) - netdata.gamma45*V(5)*cos(angleV(4) - angleV(5));
Jloss(4,4)     =     2*netdata.b45*V(5) - netdata.gamma45*V(4)*cos(angleV(4) - angleV(5));
Jloss(4,n+2)   =     netdata.gamma45*V(5)*V(4)*sin(angleV(4) - angleV(5));
Jloss(4,n+3)   =    -Jloss(4,n+2);

%Ploss46
Jloss(5,3)     =     2*netdata.b46*V(4) - netdata.gamma46*V(6)*cos(angleV(4) - angleV(6));
Jloss(5,5)     =     2*netdata.b46*V(6) - netdata.gamma46*V(4)*cos(angleV(4) - angleV(6));
Jloss(5,n+2)   =     netdata.gamma46*V(4)*V(6)*sin(angleV(4) - angleV(6));
Jloss(5,n+4)   =    -Jloss(5,n+2);


Ploss_12       =     netdata.b12*V(1)^2 + netdata.b12*V(2)^2 - netdata.gamma12*V(1)*V(2)*cos(angleV(1) - angleV(2));

Ploss_23       =     netdata.b23*V(2)^2 + netdata.b23*V(3)^2 - netdata.gamma23*V(2)*V(3)*cos(angleV(2) - angleV(3));

Ploss_34       =     netdata.b34*V(3)^2 + netdata.b34*V(4)^2 - netdata.gamma34*V(3)*V(4)*cos(angleV(3) - angleV(4));

Ploss_45       =     netdata.b45*V(4)^2 + netdata.b45*V(5)^2 - netdata.gamma45*V(4)*V(5)*cos(angleV(4) - angleV(5));

Ploss_46       =     netdata.b46*V(4)^2 + netdata.b46*V(6)^2 - netdata.gamma46*V(4)*V(6)*cos(angleV(4) - angleV(6));

Ploss = [Ploss_12;Ploss_23;Ploss_34;Ploss_45;Ploss_46];

end