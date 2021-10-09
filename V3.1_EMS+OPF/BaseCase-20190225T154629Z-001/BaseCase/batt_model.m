function [Gen_Pmax, Gen_Pmin, Soc] = batt_model( Soc, Pgen, gen_batt, k, MPCpar,Gen_Pmax,Gen_Pmin,Gen_Pmax0, Gen_Pmin0, nodi_ctrl_gen)

nb=length(gen_batt);
ncg = length(nodi_ctrl_gen);

% the battery model corresponds to a simple integrator

for b=1:nb
    
    for j=1:ncg
        
        if(nodi_ctrl_gen(j) == gen_batt(b))
          Soc(b) = Soc(b) - k(b)*Pgen(j);
        end
    end
end

Socmin = MPCpar.SocMin;
Socmax = MPCpar.SocMax;


% I compute the maximum power that is possible to ask batteries during next instant so that the
% Soc does not exceed its limits

Pmax_soc = (Soc - Socmin)./k';
Pmin_soc = (Soc - Socmax)./k';


% At this stage I check if the maximum powers imposed by the remaing SOC
% are lower with respect to the physical ones that the devices can
% nominally provide

for b=1:nb
    for j=1:ncg
                if(nodi_ctrl_gen(j)==gen_batt(b))
           Gen_Pmax(j) = min (max(Pmax_soc(b),0), Gen_Pmax0(j));
           Gen_Pmin(j) = max (min(Pmin_soc(b),0), Gen_Pmin0(j));
                end
    end
end

end