function diff = objective_function(nbd,p,TF_conc_t,RNApH_conc_t,RNApA_conc_t,mut_mat,real_data,n_strain)
% this objective function is for validation or real fitting when time-dependent RNAp is
% known

% n_strain specifies the number of strains used for fitting

% assign parameters (RNAp is not independent of energies, so it's not a
% parameter to optimize)
energyi = p(1:15);
vmax = p(16:end); % there are 24 configurations where last 2 digits are not 0

% all data follow this format:
% each row: every hour
% each column: every mutant

% generate simulated data for all WT and mutated type

if n_strain > 1
    sim_data = zeros(length(TF_conc_t),n_strain);
    parfor mm = 1:n_strain
        sim_data(:,mm) = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNApH_conc_t,RNApA_conc_t,mut_mat(mm,:),vmax);
    end
end
if n_strain == 1
    sim_data = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNApH_conc_t,RNApA_conc_t,mut_mat,vmax);
end

% now sim_data should be the same dimension as real_data

diff = weighted_msd(real_data(:),sim_data(:));

end