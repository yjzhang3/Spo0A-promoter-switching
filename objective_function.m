function diff = objective_function(nbd,p,TF_conc_t,mut_mat,real_data)

% assign parameters
RNAp_conc = p(1);
energyi = p(2:end);


% all data follow this format:
% each row: every hour
% each column: every mutant


% normalize the real data
for ll = 1:length(real_data(1,:))
    real_data(:,ll) = real_data(:,ll)./max(real_data(:,ll));
end

% generate simulated data for all WT and mutated types
sim_data = zeros(size(real_data));
for mm = 1:length(mut_mat)
    sim_data(:,mm) = time_dep_TR(nbd,energyi,TF_conc_t,RNAp_conc,mut_mat(:,mm));
end

% now sim_data should be the same dimension as real_data

diff = sum((real_data(:)-sim_data(:)).^2);

end