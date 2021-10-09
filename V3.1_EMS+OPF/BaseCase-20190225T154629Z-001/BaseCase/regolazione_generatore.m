function [Ps,Qs,dPsdV,dPsdw,dQsdV,dQsdw]=regolazione_generatore(V,w,RegG)
% Based on the voltage, the frequency and the generator characteristics
% the output powers and derivates are defined

dPsdw=0;        % Inverse RESISTIVE Droop
dQsdV=0;

Pref = RegG.Gen_Pref;
Qref = RegG.Gen_Qref;

nG = RegG.nG;

% The renewable power sources are characterized by a different droop
% configuration with respect to the fully controllable power sources.
% The renewable sources are placed at node 12 and node 13.

if ( nG ~= 12 && nG ~= 13)
    
    Qmax=RegG.Qmax;
    Qmin=RegG.Qmin;
    Pmax=RegG.Pmax;
    Pmin=RegG.Pmin;


mv = -RegG.mv;
mf =  RegG.mf;
 
% The power deviation with respect to the reference value is computed
        DP=mv*(400-V);

        
% The output active power is limited by a maximum and a minimum value

        if ((DP + Pref) >= Pmax)
                Ps = Pmax;
                dPsdV =0;

        elseif ((DP + Pref) <= Pmin)

                Ps = Pmin;
                dPsdV=0;

        else
                Ps=DP+Pref;
                dPsdV=-mv;
        end


       
        DQ = mf*(w-2*pi*50);

 % The coogenerator is characterized by a minimum power factor, therefore depending by the output active power the maximum and the
 % minimum reactive power limits may decrease. The coogenerator is placed at Node 11    
 
 if (nG == 11)      
     
     if (Qmax >= RegG.b*Ps)
                Qmax = RegG.b*Ps;
     end
     
     if (Qmin <= -RegG.b*Ps)
                Qmin = -RegG.b*Ps;
     end
 end

% The output reactive power is limited by a maximum and a minimum value

      if ((DQ + Qref) >= Qmax)
                Qs = Qmax;
                dQsdw = 0;

        elseif ((DQ + Qref) <= Qmin)

                Qs = Qmin;
                dQsdw = 0;

      else
               Qs= DQ + Qref;
               dQsdw = mf;
     end
    
else

% Renewable Sources Droop 

    Qmax=RegG.Qmax;
    Qmin=RegG.Qmin;
    
    mvs = abs(Pref)/(RegG.V(2) -RegG.V(1)); 
    
    if ((V < RegG.V(1)))
            Ps = Pref;
            dPsdV = -mvs;
            
    elseif(V >= RegG.V(1) && V < RegG.V(2))
        
            Ps = Pref - mvs*( V - RegG.V(1) );
            dPsdV = -mvs;
            
       elseif (V >= RegG.V(2))

            Ps     =     0;
            dPsdV  =     0;
            
            
    end
    
    mfs = (Qmax-Qmin)/(RegG.w(4) - RegG.w(1));
      DQ = mfs*(w-2*pi*50);

      if ((DQ + Qref) >= Qmax)
                Qs = Qmax;
                dQsdw = 0;

        elseif ((DQ + Qref) <= Qmin)

                Qs = Qmin;
                dQsdw = 0;

      else
               Qs= DQ + Qref;
               dQsdw = mfs;
     end
  
end
        
            
    

end
    


