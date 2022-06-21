function [bounds,vmax_arrays,vmax_arrayv,energyi_s,energyi_v,diff] = estimate_energy_together(inds,indv,n_strains,n_strainv,group_arrays,group_arrayv,vmax_arrays,vmax_arrayv,filename)
%% this file finds the energy for smaller set of data for *both* promoters
%% load data and choose strains
% ii tells us Pv or Ps
% n_strain tell us the number of strains to fit with
% ind tells us which strain to fit

load('promoter_activity_single.mat')

load('H_RNAP.mat', 'H_conc_t')
load('A_RNAP.mat','A_conc_t')

real_data_Ps = Ps_promoter_activity_mean(:,2:9); % exclude the real WT, instead WT is Ps only, no mutations
real_data_Ps_std = Ps_promoter_activity_std(:,2:9);
real_data_Pv = Pv_promoter_activity_mean(:,2:9);
real_data_Pv_std = Pv_promoter_activity_std(:,2:9);

TF_conc_t = get_aps(2:10);

%% assign parameters
nbd = 4;
mut_mat = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]]; 

%% set up bounds for energy

nvars = 14+numel(fieldnames(group_arrays))-numel(fieldnames(vmax_arrays))+numel(fieldnames(group_arrayv))-numel(fieldnames(vmax_arrayv));

lb_overall = zeros(nvars,1)-20;
ub_overall = zeros(nvars,1)+20;

% ub_overall([1,3]) = 5; % 0A 1 and 3 have greater affinity than 0A 2

lb_overall(4:5) = 0;
ub_overall(4:5) = 0; % keep promoter energy fixed because we know it's embedded in concentration terms
% 
% lb_overall(7) = 0;
% ub_overall(7) = 0; % no interaction energy between 0A1 and 0A3

% vmax should be positive
lb_overall(14:end) = 0;
ub_overall(14:end) = 1000;

bounds.lb = lb_overall;
bounds.ub = ub_overall;

%% index smaller dataset
    
mut_mat_new_s = mut_mat(inds,:); 
mut_mat_new_v = mut_mat(indv,:); 

real_data_new_s = real_data_Ps(:,inds);
real_data_new_v = real_data_Pv(:,indv);
%% start to fit data

[pars,diff] = fit_data_together(nbd,TF_conc_t,H_conc_t,A_conc_t,mut_mat_new_s,mut_mat_new_v,real_data_new_s,real_data_new_v,lb_overall,ub_overall,n_strains,n_strainv,...
group_arrays,group_arrayv,vmax_arrays,vmax_arrayv,filename);

%%
energyi_s = pars([1:4,6,7,9,8,10,11]);
energyi_v = pars([1:3,5,6,7,12,8,13,14]);

fng_s = fieldnames(group_arrays);
fng_v = fieldnames(group_arrayv);

fnv_s = fieldnames(vmax_arrays);
fnv_v = fieldnames(vmax_arrayv);

ol_s = numel(fnv_s); % original length of vmax array field
ol_v = numel(fnv_v); % original length of vmax array field

vmax_s = pars(15:14+numel(fng_s)-numel(fnv_s));
vmax_v = pars(15+numel(fng_s)-numel(fnv_s):end);

for k=1:numel(fng_s)
    curr_key = (fng_s{k});
    if ~isfield(vmax_arrays,curr_key)  % only if this group doesn't exist, then we have unknown vmax
        curr_value = vmax_s(k-ol_s);
        vmax_arrays.(fng_s{k}) = curr_value; % assign unknown parameters
    end
    % of the vmax's to be optimized the same as the order of the group
end

clear k

for k=1:numel(fng_v)
    curr_key = (fng_v{k});
    if ~isfield(vmax_arrayv,curr_key)  % only if this group doesn't exist, then we have unknown vmax
        curr_value = vmax_v(k-ol_v);
        vmax_arrayv.(fng_v{k}) = curr_value; % assign unknown parameters
    end
    % of the vmax's to be optimized the same as the order of the group
end

%% Ps promoter
title_name = {'Ps-WT','1*','2*','3*','12*','13*','23*','123*'};

figure();
for kk = 1:length(inds)
    subplot(4,2,kk)
    
    errorbar(TF_conc_t,real_data_Ps(:,inds(kk)),real_data_Ps_std(:,inds(kk)),'LineWidth',2,'LineStyle','none')
    hold on
    
    TR = time_dep_TR_new_wSigma(nbd,energyi_s,TF_conc_t,H_conc_t,mut_mat(inds(kk),:),vmax_arrays,group_arrays);
    plot(TF_conc_t,TR,'LineWidth',2)
    ylim([0 1000])
    xlabel('TF concentration')
    ylabel('transcription rate')
    title(string(title_name(inds(kk))))
end
%% Pv promoter
clear kk
title_name = {'Pv-WT','1*4*','2*4*','3*4*','124*','134*','234*','1234*'};
figure();
for kk = 1:length(indv)
    subplot(4,2,kk)
    
    errorbar(TF_conc_t,real_data_Pv(:,indv(kk)),real_data_Pv_std(:,indv(kk)),'LineWidth',2,'LineStyle','none')
    hold on
    
    TR = time_dep_TR_new_wSigma(nbd,energyi_v,TF_conc_t,A_conc_t,mut_mat(indv(kk),:),vmax_arrayv,group_arrayv);
    plot(TF_conc_t,TR,'LineWidth',2)
    ylim([0 1000])
    xlabel('TF concentration')
    ylabel('transcription rate')
    title(string(title_name(indv(kk))))
end
saveas(gcf,append(file(1:end-4),'.jpeg'))