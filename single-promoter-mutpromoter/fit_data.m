function [pars,M] = fit_data(nbd,TF_conc_t,RNAp_conc_t,mut_mat,real_data,lb,ub,group_vmax_array,group_delta_array,ind,file)

% input: real data, time dependent TF concentration, number of binding
% sites, time-dependent RNAP concentration, mutation matrix

% can also specify lower and upper bound

% Energyi and vmax is the parameter to find for

% can specify how many strains are grouped together to optimize

%% parameters
rng default
iter = 8;

% options = optimoptions('particleswarm','SwarmSize',50,'HybridFcn',@fmincon);

fun = @(p) objective_function(nbd,p,TF_conc_t,RNAp_conc_t,mut_mat,real_data,group_vmax_array,group_delta_array,ind);

%% particle swarm
nvars = 10+numel(fieldnames(group_vmax_array))+numel(fieldnames(group_delta_array));
pars_all = zeros(iter,nvars);
diff_all = zeros(iter,1);
for n = 1:iter
    [x,fval,~,~] = particleswarm(fun,nvars,lb,ub);
    pars_all(n,:) = x;
    diff_all(n) = fval;
    create_parameter_file(file, x, fval);
end

[M,I] = min(diff_all);
pars =  pars_all(I,:);
end