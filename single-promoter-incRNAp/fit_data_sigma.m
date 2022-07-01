function [pars,M] = fit_data_sigma(nbd,energyi,TF_conc_t,mut_RNAp,real_data,vmax_array,group_array,file)

% input: real data, time dependent TF concentration, number of binding
% sites, mutation type 

% will not specify lower and upper bound because it's hard to know

% time-dependent [RNAp] is the parameter to find for

% we can just let all energies being 1 and let [RNAP] cinlude the promoter
% energy. 

%% parameters
rng default
iter = 8;

fnv = fieldnames(vmax_array);

if isempty(fnv)
    nvars = length(TF_conc_t)+numel(fieldnames(group_array));
end
if ~isempty(fnv)
    nvars = length(TF_conc_t)+numel(fieldnames(group_array))-numel(fieldnames(vmax_array));
end

lb = zeros(nvars,1);
ub = zeros(nvars,1)+25;
ub(length(TF_conc_t)+1:end) = 1000; % upper bound for vmax is 1000

fun = @(p) objective_function_sigma(nbd,energyi,TF_conc_t,p,mut_RNAp,real_data,vmax_array,group_array);


%% particle swarm
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