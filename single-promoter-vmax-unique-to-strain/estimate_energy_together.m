function [pars,diff] = estimate_energy_together(group_array_s,group_array_v,ind_s,ind_v,filename)
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

nvars = 14+numel(fieldnames(group_array_s))+numel(fieldnames(group_array_v));

lb_overall = zeros(nvars,1)-30;
ub_overall = zeros(nvars,1)+30;

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
    
mut_mat_new_s = mut_mat(ind_s,:); 
mut_mat_new_v = mut_mat(ind_v,:); 

real_data_new_s = real_data_Ps(:,ind_s);
real_data_new_v = real_data_Pv(:,ind_v);
%% start to fit data

[pars,diff] = fit_data_together(nbd,TF_conc_t,H_conc_t,A_conc_t,mut_mat_new_s,mut_mat_new_v,real_data_new_s,real_data_new_v,lb,ub,...
    group_array_s,group_array_v,ind_s,ind_v,filename);

%% assign group array and vmax array
[~,energyi_s,energyi_v,vmax_per_strain_final_s,vmax_per_strain_final_v] = ...
    objective_function_together(nbd,pars,TF_conc_t,H_conc_t,A_conc_t,...
    mut_mat_new_s,mut_mat_new_v,real_data_new_s,real_data_new_v,...
    group_array_s,group_array_v,ind_s,ind_v);

%% Ps promoter
title_name = {'Ps-WT','1*','2*','3*','12*','13*','23*','123*'};

figure();
for kk = 1:length(ind_s)
    subplot(4,2,kk)
    
    errorbar(TF_conc_t,real_data_Ps(:,ind_s(kk)),real_data_Ps_std(:,ind_s(kk)),'LineWidth',2,'LineStyle','none')
    hold on
    
    TR = time_dep_TR_new_wSigma(nbd,energyi_s,TF_conc_t,H_conc_t,mut_mat(ind_s(kk),:),vmax_per_strain_final_s(kk));
    plot(TF_conc_t,TR,'LineWidth',2)
    ylim([0 1000])
    xlabel('TF concentration')
    ylabel('transcription rate')
    title(string(title_name(ind_s(kk))))
end
%% Pv promoter
clear kk
title_name = {'Pv-WT','1*4*','2*4*','3*4*','124*','134*','234*','1234*'};
figure();
for kk = 1:length(ind_v)
    subplot(4,2,kk)
    
    errorbar(TF_conc_t,real_data_Pv(:,ind_v(kk)),real_data_Pv_std(:,ind_v(kk)),'LineWidth',2,'LineStyle','none')
    hold on
    
    TR = time_dep_TR_new_wSigma(nbd,energyi_v,TF_conc_t,A_conc_t,mut_mat(ind_v(kk),:),vmax_per_strain_final_v(kk));
    plot(TF_conc_t,TR,'LineWidth',2)
    ylim([0 1000])
    xlabel('TF concentration')
    ylabel('transcription rate')
    title(string(title_name(ind_v(kk))))
end