function diff = objective_function_sigma(nbd,energyi,TF_conc_t,p,mut_RNAp,real_data)
% this objective function is for estimating how RNAp (sigma factor
% specifically) changes with time

% assign parameters 
RNAp_conc_t = p(1:end-1);
vmax = p(end);

% all data follow this format:
% each row: every hour
% each column: every mutant

% generate simulated data for all WT and mutated types
sim_data = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNAp_conc_t,mut_RNAp,vmax);

% now sim_data should be the same dimension as real_data

diff = weighted_msd(real_data(:),sim_data(:));

end