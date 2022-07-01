function [pars,M] = fit_data(nbd,TF_conc_t,RNAp_conc_t,mut_mat,real_data,lb,ub,n_strain,group_array,vmax_array,file)
% input: real data, time dependent TF concentration, number of binding
% sites, time-dependent RNAP concentration, mutation matrix

% can also specify lower and upper bound

% Energyi and vmax is the parameter to find for

% can specify how many strains are grouped together to optimize

%% parameters
rng default
iter = 8;

% options = optimoptions('particleswarm','SwarmSize',50,'HybridFcn',@fmincon);

fun = @(p) objective_function(nbd,p,TF_conc_t,RNAp_conc_t,mut_mat,real_data,n_strain,vmax_array,group_array);

fnv = fieldnames(vmax_array);

if isempty(fnv)
    nvars = 10+numel(fieldnames(group_array));
end
if ~isempty(fnv)
    nvars = 10+numel(fieldnames(group_array))-numel(fieldnames(vmax_array));
end

% number of unknown vmax is the difference in the number of fields between
% group_array and vmax_array

% energy of each binding site and vmax of each configuration

%% particle swarm
pars_all = zeros(iter,nvars);
diff_all = zeros(iter,1);
for n = 1:iter
    [x,fval,~,~] = particleswarm(fun,nvars,lb,ub);
%     [x,fval,~,~] = fmincon(fun,zeros(nvars,1),[],[],[],[],lb,ub);
%     [x,fval,~,~] = fminsearch(fun,lb);
    pars_all(n,:) = x;
    diff_all(n) = fval;
    create_parameter_file(file, x, fval);
end

[M,I] = min(diff_all);
pars =  pars_all(I,:);
end