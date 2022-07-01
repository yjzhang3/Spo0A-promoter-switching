function diff = objective_function_sigma(nbd,energyi,TF_conc_t,p,mut_RNAp,real_data,vmax_array,group_array)
% this objective function is for estimating how RNAp (sigma factor
% specifically) changes with time

% assign parameters 
RNAp_conc_t = p(1:length(TF_conc_t));
vmax_raw = p(length(TF_conc_t)+1:end);

%%

% create the vmax_array here:
% I would always put known vmax as g1!!! So this can avoid confusion and
% don't worry about assigning order
fn = fieldnames(group_array);
fnv = fieldnames(vmax_array);

if isempty(fnv)
    ol = 0;
end
if ~isempty(fnv)
    ol = numel(fieldnames(vmax_array)); % original length of vmax array field
end


for k=1:numel(fn)
    curr_key = (fn{k});
    if ~isfield(vmax_array,curr_key)  % only if this group doesn't exist, then we have unknown vmax
        curr_value = vmax_raw(k-ol);
        vmax_array.(fn{k}) = curr_value; % assign unknown parameters to the corresponding vmax group
    end
    % of the vmax's to be optimized the same as the order of the group
end

%%
% all data follow this format:
% each row: every hour
% each column: every mutant

% generate simulated data for all WT and mutated type
n_strain = length(mut_RNAp(:,1));

if n_strain > 1
    sim_data = zeros(length(TF_conc_t),n_strain);
    for mm = 1:n_strain
        sim_data(:,mm) = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNAp_conc_t,mut_RNAp(mm,:),vmax_array,group_array);
    end
end
if n_strain == 1
    sim_data = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNAp_conc_t,mut_RNAp,vmax_array,group_array);
end



% now sim_data should be the same dimension as real_data

diff = weighted_msd(real_data(:),sim_data(:));

end