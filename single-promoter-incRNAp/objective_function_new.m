function diff = objective_function_new(nbd,p,TF_conc_t,RNAp_conc_t,mut_mat,real_data)
% this objective function is for validation or real fitting when time-dependent RNAp is
% known

% assign parameters (RNAp is not independent of energies, so it's not a
% parameter to optimize)
energyi = p(1:end-1);
vmax = p(end);

% all data follow this format:
% each row: every hour
% each column: every mutant

% generate simulated data for all WT and mutated types
sim_data = zeros(length(mut_mat),length(TF_conc_t));
for mm = 1:length(mut_mat)
    sim_data(mm,:) = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNAp_conc_t,mut_mat(mm,:),vmax);
end

% now sim_data should be the same dimension as real_data

diff = weighted_msd(real_data(:),sim_data(:));

end