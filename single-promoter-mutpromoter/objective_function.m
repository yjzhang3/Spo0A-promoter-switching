function [diff,sim_data] = objective_function(nbd,p,TF_conc_t,RNAp_conc_t,mut_mat,real_data,group_vmax_array,group_delta_array,ind)

% this objective function is for validation or real fitting when time-dependent RNAp is
% known
n_strain = length(ind);
% n_strain specifies the number of strains used for fitting

% assign parameters (RNAp is not independent of energies, so it's not a
% parameter to optimize)

energyi = p(1:10);
delta_raw = p(11:10+numel(fieldnames(group_delta_array)));
vmax_raw = p(11+numel(fieldnames(group_delta_array)):end);

fnv = fieldnames(group_vmax_array);
vmax_per_strain = zeros(8,1);
for k=1:numel(fnv)
    curr_ind = group_vmax_array.(fnv{k});
    curr_vmax = vmax_raw(k);
    vmax_per_strain(curr_ind) = curr_vmax;
end
vmax_per_strain_final = vmax_per_strain(ind);

clear k
energyi_per_strain = repmat(energyi, 8, 1); % potentially 8 different set of energies
fnd = fieldnames(group_delta_array);
for k=1:numel(fnd)
    curr_ind = group_delta_array.(fnd{k});
    curr_delta = delta_raw(k);
    energyi_per_strain(curr_ind,4) = curr_delta; % change the promoter affinity to a different parameter
end
energyi_per_strain_final = energyi_per_strain(ind,:);

% all data follow this format:
% each row: every hour
% each column: every mutant

% generate simulated data for all WT and mutated type

if n_strain > 1
    sim_data = zeros(length(TF_conc_t),n_strain);
    for mm = 1:n_strain
        sim_data(:,mm) = time_dep_TR_new_wSigma(nbd,energyi_per_strain_final(mm,:),TF_conc_t,RNAp_conc_t,mut_mat(mm,:),vmax_per_strain_final(mm));
    end
end
% if n_strain == 1
%     sim_data = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNAp_conc_t,mut_mat,vmax_raw);
% end

% now sim_data should be the same dimension as real_data

diff = weighted_msd(real_data(:),sim_data(:));

end