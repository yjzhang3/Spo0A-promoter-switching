function p = prob_transcription(i,Ni,Nr,energy,concSpo0A,concRNAP)
% output: probability of each state 
% input: an array of configurations available after one or more binding site is blocked,
% binding affinity e for all configurations, as well as the concentration of Spo0A 
% dimer and RNApolymerase 

% constant in Bolzmann distribution
kT = 1;

% number of Spo0A dimers and RNA bound for each configuration
Ni = [1,1,1,1,2,2,2,3];
Nr = [0,1,1,1,1,1,1,1];

% denominator (partition function Z)
deno = partition_fun(energy,concSpo0A,concRNAP);

% numerator
nume = 0;
for ll=1:length(Ni)
    if ll ~= i % exclude the configuration that is bound
        nume = nume+exp(-energy(ll)/kT)*concSpo0A^Ni(ll)*concRNAP^Nr(ll);
    end 
end

p = nume/deno;

end




