function TR = transcription_rate(nbd,energyi,TF_conc,RNAp_conc)
% input: 
% nbd: # of binding sites
% TF_conc and mRNA conc (assuming only one TF involved in the problem)
% energyi: binding energy of each binding site and promoter


% output: transcription rate

% generate time dependent concentration


% partition function Z_on
Zon = Z_on(nbd,energyi,TF_conc,RNAp_conc);

% partition function Z_off
Zoff = Z_off(nbd,energyi,TF_conc,RNAp_conc);
    
TR = Zon/(Zon+Zoff);

end