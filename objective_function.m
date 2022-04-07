function diff = objective_function(nbd,p,TF_conc_t,RNAp_conc,mut_mat,real_data)

% assign parameters
energyi = p;


% all data follow this format:
% each row: every hour
% each column: every mutant


% normalize the real data
for ll = 1:length(real_data(1,:))
    real_data(:,ll) = real_data(:,ll)./max(real_data(:,ll));
end

% generate simulated data for all WT and mutated types
sim_data = zeros(length(mut_mat),length(TF_conc_t));
parfor mm = 1:length(mut_mat)
    sim_data(:,mm) = time_dep_TR_new(nbd,energyi,TF_conc_t,RNAp_conc,mut_mat(mm,:));
end

% now sim_data should be the same dimension as real_data

diff = weighted_msd(real_data(:),sim_data(:));

end