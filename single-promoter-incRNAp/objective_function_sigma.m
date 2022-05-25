function diff = objective_function_sigma(nbd,energyi,TF_conc_t,p,mut_RNAp,real_data)
% this objective function is for estimating how RNAp (sigma factor
% specifically) changes with time

% assign parameters 
RNAp_conc_t = p(1:length(TF_conc_t));
vmax_raw = p(end);


% create the vmax_array here:
% this is a special case for fitting 123* data alone
vmax_array.g1 = 0;
vmax_array.g2 = vmax_raw;

group_array.g1 = [2:8]; % with respect to the on_config created by the algorithm,
% 2:8 are strains with some mutations but not all
group_array.g2 = 1;

% all data follow this format:
% each row: every hour
% each column: every mutant

% generate simulated data for all WT and mutated types
sim_data = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNAp_conc_t,mut_RNAp,vmax_array,group_array);

% now sim_data should be the same dimension as real_data

diff = weighted_msd(real_data(:),sim_data(:));

end