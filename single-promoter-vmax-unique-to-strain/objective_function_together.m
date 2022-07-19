function [diff,energyi_s,energyi_v,vmax_per_strain_final_s,vmax_per_strain_final_v] = ...
    objective_function_together(nbd,p,TF_conc_t,H_conc_t,A_conc_t,...
    mut_mat_s,mut_mat_v,real_data_s,real_data_v,group_array_s,group_array_v,...
    ind_s,ind_v)
% this objective function is for validation or real fitting when time-dependent RNAp is
% known

% n_strain specifies the number of strains used for fitting
n_strain_s = length(ind_s);
n_strain_v = length(ind_v);

% assign parameters (RNAp is not independent of energies, so it's not a
% parameter to optimize)
% 1,2,3,s,v,12,13,23,1s,2s,3s,1v,2v,3v
single_energyi = p(1:5);
TF_int_energyi = p(6:8);
Ps_int_energyi = p(9:11);
Pv_int_energyi = p(12:14);

% 1,2,3,p,12,13,1p,23,2p,3p
energyi_s = zeros(10,1);
energyi_s(1:4) = single_energyi(1:4);
energyi_s([5,6,8]) = TF_int_energyi;
energyi_s([7,9,10]) = Ps_int_energyi;

energyi_v = zeros(10,1);
energyi_v(1:4) = single_energyi([1:3,5]);
energyi_v([5,6,8]) = TF_int_energyi;
energyi_v([7,9,10]) = Pv_int_energyi;

%% create vmax array
fng_s = fieldnames(group_array_s);
fng_v = fieldnames(group_array_v);
vmax_s = p(15:14+numel(fng_s)-numel(fnv_s));
vmax_v = p(15+numel(fng_s)-numel(fnv_s):end);

vmax_per_strain_s = zeros(8,1);
vmax_per_strain_v = zeros(8,1);
for k=1:numel(fng_s)
    curr_ind = group_array_s.(fng_s{k});
    curr_vmax = vmax_s(k);
    vmax_per_strain_s(curr_ind) = curr_vmax;
end
vmax_per_strain_final_s = vmax_per_strain_s(ind_s);

clear k

for k=1:numel(fng_v)
    curr_ind = group_array_v.(fng_v{k});
    curr_vmax = vmax_v(k);
    vmax_per_strain_v(curr_ind) = curr_vmax;
end
vmax_per_strain_final_v = vmax_per_strain_s(ind_v);
%% 

% all data follow this format:
% each row: every hour
% each column: every mutant

% generate simulated data for all WT and mutated type

sim_data_s = zeros(length(TF_conc_t),n_strains);
sim_data_v = zeros(length(TF_conc_t),n_strainv);

for mm = 1:n_strain_s
    sim_data_s(:,mm) = time_dep_TR_new_wSigma(nbd,energyi_s,TF_conc_t,H_conc_t,mut_mat_s(mm,:),vmax_per_strain_final_s(mm));
end

for nn = 1:n_strain_v
    sim_data_v(:,nn) = time_dep_TR_new_wSigma(nbd,energyi_v,TF_conc_t,A_conc_t,mut_mat_v(nn,:),vmax_per_strain_final_v(mm));
end

sim_data = [sim_data_s;sim_data_v];
real_data = [real_data_s;real_data_v];

% now sim_data should be the same dimension as real_data

diff = weighted_msd(real_data(:),sim_data(:));

end