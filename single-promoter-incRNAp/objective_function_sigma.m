function diff = objective_function_sigma(nbd,p,TF_conc_t,mut_RNAp,real_data)
% this objective function is for estimating how RNAp (sigma factor
% specifically) changes with time

% assign parameters 
energy_len = nbd+nbd*(nbd-1)/2;
energyi = p(1:energy_len);
RNAp_conc_t = p(energy_len+1:end);

% all data follow this format:
% each row: every hour
% each column: every mutant

% generate simulated data for all WT and mutated types
sim_data = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNAp_conc_t,mut_RNAp);

% now sim_data should be the same dimension as real_data

diff = weighted_msd(real_data(:),sim_data(:));

end