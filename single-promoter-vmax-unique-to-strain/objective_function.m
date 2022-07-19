function [diff,vmax_per_strain_final] = objective_function(nbd,p,TF_conc_t,RNAp_conc_t,mut_mat,real_data,group_array,ind)

% this objective function is for validation or real fitting when time-dependent RNAp is
% known
n_strain = length(ind);
% n_strain specifies the number of strains used for fitting

% assign parameters (RNAp is not independent of energies, so it's not a
% parameter to optimize)

energyi = p(1:10);
vmax_raw = p(11:end);

fn = fieldnames(group_array);
vmax_per_strain = zeros(8,1);
for k=1:numel(fn)
    curr_ind = group_array.(fn{k});
    curr_vmax = vmax_raw(k);
    vmax_per_strain(curr_ind) = curr_vmax;
end
vmax_per_strain_final = vmax_per_strain(ind);

% all data follow this format:
% each row: every hour
% each column: every mutant

% generate simulated data for all WT and mutated type

if n_strain > 1
    sim_data = zeros(length(TF_conc_t),n_strain);
    for mm = 1:n_strain
        sim_data(:,mm) = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNAp_conc_t,mut_mat(mm,:),vmax_per_strain_final(mm));
    end
end

diff = weighted_msd(real_data(:),sim_data(:));

end