% Power Flow di reti elettriche connesse ad un nodo di saldo (rete
% prevalente), partendo dai parametri di rete, generatori e carichi
% 
%   results=connected_PF(netdata,RegG,RegL)
% 
% La funzione restituisce:
% - results.V  - Vettori spaziali (Park) delle tensioni nodali [V]
% - results.I  - Vettori spaziali (Park) delle correnti di linea [A]
% - results.S  - Potenze complesse in ingresso nei nodi di rete [VA]
% 
% Gli input sono:
% - dati della rete
%     netdata.ni          - nodi iniziali delle linee
%     netdata.nf          - nodi finali delle linee
%     netdata.ns          - nodo di saldo/riferimento della rete
%     netdata.Vs          - tensione del nodo di riferimento [V]
%     netdata.wn          - frequenza a cui opera della rete [rad/s]
%     netdata.R/L/C/G/T/k - resistenze/induttanze/capacità/conduttanze/
%                           elastanze/rapporti di trasformazione delle
%                           linee [Ohm / H / F / 1/Ohm / 1/H / V/V]
% 
% - Parametri di regolazione dei generatori
%     RegG...             - parametri elaborati dalla funzione
%                           'regolazione_generatore'
%   La funzione di regolazione dei generatori deve essere nella forma:
%      [P,Q,dPdV,dPdw,dQdV,dQdw]=regolazione_generatore(V,w,RegG)
%      - V     è il modulo della tensione ai morsetti del generatore
%      - w     è la frequenza di rete
%      - P     è la potenza attiva erogata dal generatore
%      - Q     è la potenze reattiva scambiata dal generatore (positiva se
%              in sovraeccitazione)
%      - dPdV  è la derivata della potenza attiva rispetto alla tensione
%      - dPdw  è la derivata della potenza attiva rispetto alla frequenza
%      - dQdV  è la derivata della potenza reattiva rispetto alla tensione
%      - dQdw  è la derivata della potenza reattiva rispetto alla frequenza
% 
% - Parametri di regolazione dei carichi
%     RegL...             - parametri elaborati dalla funzione
%                           'regolazione_carico'
%   La funzione di regolazione dei carichi deve essere nella forma:
%      [P,Q,dPdV,dPdw,dQdV,dQdw]=regolazione_carico(V,w,RegG)
%      - V     è il modulo della tensione ai morsetti del generatore
%      - w     è la frequenza di rete
%      - P     è la potenza attiva consumata dal generatore
%      - Q     è la potenze reattiva scambiata dal carico (positiva se di
%              natura induttiva)
%      - dPdV  è la derivata della potenza attiva rispetto alla tensione
%      - dPdw  è la derivata della potenza attiva rispetto alla frequenza
%      - dQdV  è la derivata della potenza reattiva rispetto alla tensione
%      - dQdw  è la derivata della potenza reattiva rispetto alla frequenza
% 
% Opzioni aggiuntive possono essere trasmesse/estratte alla/dalla funzione
% utilizzando la seguente sintassi
% 
%   [results,details]=connected_PF(netdata,RegG,RegL,options)
% 
% dove, facoltativamente, possono essere indicati i seguenti valori:
% - options.maxD       - Massime variazioni (per iterazione) dei valori
%                        di modulo tensioni, fase tensioni e frequenza.
%                        Di default, maxD=[0.20,0.10], a cui corrispondono
%                           20%       variazione massima modulo di tensione
%                           0.1 rad   variazione massima fase tensione
% - options.IterMax    - Massimo numero di iterazioni di Power Flow per
%                        per raggiungere la convergenza (di default
%                        IterMax=50);
% - options.IterOK     - Numero di iterazioni di Power Flow all'interno
%                        della tolleranza eseguite le quali il sistema si
%                        considera in convergenza (di default IterOK=2)
% - options.Tol        - Tolleranze per cui si considera raggiunta la
%                        convergenza. Di default si ha Tol=[1e-5,0.1] per
%                        cui il sistema si ritiene a convergenza quando la
%                        variazione di tensione/frequenza/potenza tra due
%                        iterazioni successive è inferiore a
%                           0.001%    variazione modulo di tensione
%                           0.10 VA   variazione potenza
% - options.Vin        - Valori iniziali delle tensioni nodali, per cui
%                        viene eseguita la prima iterazione di Power Flow.
%                        Di default vengono calcolate le tensioni della
%                        rete a vuoto [V]
% 
% La funzione restituisce dettagli riguardo al processo di Power flow. In
% particolare:
% - details.dv/ds        - valori dei residui di tensione/potenza per ogni
%                          step del Power Flow
% - details.iterations   - numero di iterazioni effettuate per trovare la
%                          soluzione
% - details.convergence  - indica la raggiunta convergenza (1 = converge)
% - details.Yn           - matrice delle ammettenze nodali [1/Ohm]
% - details.Sload        - Potenze complesse assorbite dai carichi [VA]
% - details.Sgen         - Potenze complesse iniettate dai generatori [VA]
