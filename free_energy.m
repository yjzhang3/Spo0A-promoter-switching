function Gtot = free_energy(conc_arr,num_arr)
% input: 
% conc_arr: concentration of each protein present i in this configuration
% num_arr: number of each protein i present in this config
% same index refers to the same type of protein
% right now, assume only one type of protein can bind to one binding box,
% so that these arrays should be the same length as the number of binding
% sites. Each position in the array represent the concentration of TF and
% number TF bound at this particular binding box. 

% output: 
% Gtot: free energy of each *configuration*

%% constant
% standard free energy of the configuration
G0 = 20;

%% formula
Gtot = G0 + sum(num_arr.*log(conc_arr+0.000000001));

end