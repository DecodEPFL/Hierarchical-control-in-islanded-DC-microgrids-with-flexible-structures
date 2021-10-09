function [Ps,Qs,dPsdV,dPsdw,dQsdV,dQsdw]=regolazione_generatore(V,w,RegG)
% calcola la potenza attiva e la potenza reattiva per i valori di tensione
% e frequenza

Ps = RegG.P;
Qs = RegG.Q;
dPsdV=0;
dPsdw=0;
dQsdV=0;
dQsdw=0;
end
    


