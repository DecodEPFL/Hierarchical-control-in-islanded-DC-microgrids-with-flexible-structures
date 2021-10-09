function [Ps,Qs,dPsdV,dPsdw,dQsdV,dQsdw]=regolazione_carico(V,w,RegL)
% calcola la potenza attiva e la potenza reattiva per i valori di tensione
% e frequenza

Ps = RegL.P;
Qs = RegL.Q;
dPsdV=0;
dPsdw=0;
dQsdV=0;
dQsdw=0;
end
