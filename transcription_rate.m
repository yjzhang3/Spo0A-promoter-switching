function TR = transcription_rate(nbd,TF_conc,energyi,typei,coopi)
% input: 
% nbd: # of binding sites
% TF_concentration (assuming only one TF involved in the problem)
% energyi: binding energy of each binding site and promoter
% typei: active or repressive power of each TF (>1 or <1, TF-mRNA
% interaction)
% coopi: cooperativity of each configuration (should have 2^(nbd-1) terms
% because we only consider cooperativity between each TF

% output: transcription rate

% generate time dependent concentration


% partition function Z_on
Zon = Z_on(nbd,energyi,typei,coopi,TF_conc);

% partition function Z_off
Zoff = Z_off(nbd,energyi,coopi,TF_conc);
    
TR = Zon/(Zon+Zoff);

end