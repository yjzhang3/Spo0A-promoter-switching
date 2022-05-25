function diff = objective_function(nbd,p,TF_conc_t,RNAp_conc_t,mut_mat,real_data,n_strain,group_array)
% this objective function is for validation or real fitting when time-dependent RNAp is
% known

% n_strain specifies the number of strains used for fitting

% assign parameters (RNAp is not independent of energies, so it's not a
% parameter to optimize)
energyi = p(1:10);
vmax_raw = p(11:end);

% create the vmax_array here:
fn = fieldnames(group_array);
for k=1:numel(fn)
    vmax_array.(fn{k}) = vmax_raw(k); % this is only possible if the order
    % of the vmax's to be optimized the same as the order of the group
end

% all data follow this format:
% each row: every hour
% each column: every mutant

% generate simulated data for all WT and mutated type

if n_strain > 1
    sim_data = zeros(length(TF_conc_t),n_strain);
    parfor mm = 1:n_strain
        sim_data(:,mm) = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNAp_conc_t,mut_mat(mm,:),vmax_array,group_array);
    end
end
if n_strain == 1
    sim_data = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNAp_conc_t,mut_mat,vmax_array,group_array);
end

% now sim_data should be the same dimension as real_data

diff = weighted_msd(real_data(:),sim_data(:));

end