function [Ps,Qs,dPsdV,dPsdw,dQsdV,dQsdw]=regolazione_carico(V,w,RegL)
% Based on the voltage, the frequency and on the load characteristic the
% absorbed active and reactive power and the relative derivatives are
% computed

%Nominal value of voltages and frequency
Vn=400;
wn=2*pi*50;


nL=RegL.nL;

% Loads at Node 2 and Node 6 are Parallel RLC Loads

if ((nL == 2) || (nL == 6))

    Pr   = RegL.Pr;
    Qind = RegL.Qind;
    Qcap = RegL.Qcap;
    
    R =  Vn*Vn/Pr;
    L = Vn*Vn/(wn*Qind);
    C =  -Qcap/(Vn*Vn*wn);
    
    Ps = -V*V/R;
    Qs = -V*V/(w*L) + V*V*w*C;
    
    dPsdV = -2*V/R;
    dPsdw = 0;
    
    dQsdV =   - 2*V/(w*L) + 2*V*C*w;
    dQsdw =   +V*V/(L*w*w) + V*V*C;
    
% Loads at Node 5, Node 9 and Node 11 are costant loads

elseif (nL == 5 || nL==9 || nL==11)
        
        Ps=-real(RegL.Sn);
        Qs=-imag(RegL.Sn);
        
        dPsdV = 0;
        dPsdw = 0;
        dQsdV = 0;
        dQsdw = 0;
   
end
  

end


