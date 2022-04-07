function pars = fit_data_new(nbd,TF_conc_t,RNAp_conc,mut_mat,real_data,lb,ub)
% input: real data, time dependent TF concentration, number of binding
% sites, RNAP concentration, mutation matrix
% files with "_new" is for validation

% can also specify lower and upper bound

% Energyi is the parameter to find for

%% parameters
rng default
iter = 10;

fun = @(p) objective_function_new(nbd,p,TF_conc_t,RNAp_conc,mut_mat,real_data);

nvars = 10;

%% particle swarm
pars = zeros(iter,nvars);
parfor n = 1:iter
    x = particleswarm(fun,nvars,lb,ub);
    pars(n,:) = x;
end

end