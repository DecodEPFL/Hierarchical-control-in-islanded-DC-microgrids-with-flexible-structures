function volt = ZIPloads_par2(t,y,net_par)

PL  = net_par.PL;
I   = net_par.I;
X   = net_par.X;
N   = net_par.N;

I(1:N)=I(1:N)+(PL)*y(1:N).^(-1);

volt=X*y+I;

end