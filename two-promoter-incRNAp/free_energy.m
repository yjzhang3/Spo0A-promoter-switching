function Gtot = free_energy(config,energyi,TF_conc,RNAp_conc)
% input: 
% conc_arr: concentration of each protein i present in this configuration
% num_arr: number of each protein i present in this config

% same index refers to the same type of protein

% it's NOT necessarily the case that the concentration and number of a
% particular protein is unique to each binding box. So the length of these
% two arrays has nothing to do with each binding box. We simply use these
% two arrays to represent the number and concentration of a particular
% protein

% let's be consistent, the last element of both arrays represent mRNA; all
% the other elements before this refer to TF
% output: 
% Gtot: free energy of each *configuration*

%% standard free energy of the configuration
G0 = stand_energy(config,energyi);


%% concentration and number array for this configuration

[num_arr,conc_arr] = gen_config_prof(config,TF_conc,RNAp_conc);


%% formula 
Gtot = G0 - sum(num_arr.*log(conc_arr)); 

% Gtot = G0; % start with G0
% for cc = 1:length(conc_arr)
%     if conc_arr(cc) ~= 0
%         ct = num_arr(cc)*log(conc_arr(cc)); % correction term at this iteration
%     else
%         ct = 0;
%     end
%     
%     Gtot = Gtot-ct;
% end

end