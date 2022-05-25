function [pars,M] = fit_data_sigma(nbd,energyi,TF_conc_t,mut_RNAp,real_data)
% input: real data, time dependent TF concentration, number of binding
% sites, mutation type for 123*

% will not specify lower and upper bound because it's hard to know

% time-dependent [RNAp] is the parameter to find for

% we can just let all energies being 1 and let [RNAP] cinlude the promoter
% energy. 

%% parameters
rng default
iter = 4;
n_pars = length(TF_conc_t)+1; % all time-dependent RNAp plus a vmax
lb = zeros(n_pars,1);
ub = zeros(n_pars,1)+100;
ub(end) = 100000;

fun = @(p) objective_function_sigma(nbd,energyi,TF_conc_t,p,mut_RNAp,real_data);

nvars = n_pars;

%% particle swarm
pars_all = zeros(iter,nvars);
diff_all = zeros(iter,1);
for n = 1:iter
    [x,fval,~,~] = particleswarm(fun,nvars,lb,ub);
    pars_all(n,:) = x;
    diff_all(n) = fval;
end

[M,I] = min(diff_all);
pars =  pars_all(I,:);
end