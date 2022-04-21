function [pars,M] = fit_data_sigma(nbd,TF_conc_t,mut_RNAp,real_data,lb,ub)
% input: real data, time dependent TF concentration, number of binding
% sites, mutation type for 123*

% can also specify lower and upper bound

% time-dependent [RNAp] and energyi is the parameter to find for

%% parameters
rng default
iter = 5;

fun = @(p) objective_function_sigma(nbd,p,TF_conc_t,mut_RNAp,real_data);

nvars = 16;

%% particle swarm
pars_all = zeros(iter,nvars);
diff_all = zeros(iter,1);
for n = 1:iter
    x = particleswarm(fun,nvars,lb,ub);
    pars_all(n,:) = x;
    diff_all(n) = fun(x);
end

[M,I] = min(diff_all);
pars =  pars_all(I,:);
end