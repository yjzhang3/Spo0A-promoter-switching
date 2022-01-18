function Gtot = free_energy(conc_arr,num_arr)
% input: 
% conc_arr: concentration of each protein present i in this configuration
% num_arr: number of each protein i present in this config
% same index refers to the same type of protein

% output: 
% Gtot: free energy of each *configuration*

%% constant
% standard free energy of the configuration
G0 = 20;

%% formula
Gtot = G0 + sum(num_arr.*log(conc_arr));

end