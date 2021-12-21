function Z = partition_fun(energy,concSpo0A,concRNAP)
% output: the partition function of our system
% input: arry of energy values, concentration of Spo0A dimers and RNAPolymerase

% constant
kT = 1;

% number of Spo0A dimers and RNA bound for each configuration
Ni = [1,1,1,1,2,2,2,3];
Nr = [0,1,1,1,1,1,1,1];


Z = 0;
for ll = 1:length(Ni)
    Z = Z+exp(-energy(ll)/kT)*concSpo0A^Ni(ll)*concRNAP^Nr(ll);
end 