function [pars,diff] = estimate_energy_together(n_strains,n_strainv,inds,indv)
%% this file finds the energy for smaller set of data
%% load data and choose strains
% ii tells us Pv or Ps
% n_strain tell us the number of strains to fit with
% ind tells us which strain to fit

load('promoter_activity_single.mat')

load('sigmaA_conc_t.mat', 'A_conc_t')
load('sigmaH_conc_t.mat','H_conc_t')

real_data_Ps = Ps_promoter_activity_mean(:,2:9); % exclude the real WT, instead WT is Ps only, no mutations
real_data_Ps_std = Ps_promoter_activity_std(:,2:9);
real_data_Pv = Pv_promoter_activity_mean(:,2:9);
real_data_Pv_std = Pv_promoter_activity_std(:,2:9);

TF_conc_t = get_aps(2:10);

%% assign parameters
nbd = 4;
mut_mat = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[0,0,1];[0,0,0]];

%% set up bounds for energy

lb_overall = zeros(32,1)-20;
ub_overall = zeros(32,1)+20;

% lb_overall(1:3) = 1; % binding affinity of each site > 1 (nonzero and greater than promoter affinity)
% lb_overall(7) = 5; % must have one repressor

lb_overall(4) = 0;
ub_overall(4) = 0; % keep promoter energy fixed because we know it's embedded in concentration terms

% vmax should be positive
lb_overall(17:end) = 0;
ub_overall(17:end) = 10000;

%% index smaller dataset
    
mut_mat_new_s = mut_mat(inds,:); 
mut_mat_new_v = mut_mat(indv,:); 

real_data_new_s = real_data_Ps(:,inds);
real_data_new_v = real_data_Pv(:,indv);
%% start to fit data

[pars,diff] = fit_data_together(nbd,TF_conc_t,H_conc_t,A_conc_t,mut_mat_new_s,mut_mat_new_v,real_data_new_s,real_data_new_v,lb_overall,ub_overall,n_strains,n_strainv);

%%
energyi_s = pars(1:10);
energyi_v = pars([1:4,11:16]);
vmax_s = pars(17:24);
vmax_v = pars(25:32);

%% Ps promoter
title_name = {'Pv','1*4*','2*4*','3*4*','124*','134*','234*','1234*'};
figure();
for kk = 1:length(inds)
    subplot(4,2,kk)
    
    errorbar(TF_conc_t,real_data_Ps(:,inds(kk)),real_data_Ps_std(:,inds(kk)),'LineWidth',4)
    hold on
    
    TR = time_dep_TR_new_wSigma(nbd,energyi_s,TF_conc_t,H_conc_t,mut_mat(inds(kk),:),vmax_s);
    plot(TF_conc_t,TR,'+','LineWidth',2)
    xlabel('TF concentration')
    ylabel('transcription rate')
    title(string(title_name(inds(kk))))
end
%% Pv promoter
title_name = {'Ps','1*','2*','3*','12*','13*','23*','123*'};
figure();
for kk = 1:length(indv)
    subplot(4,2,kk)
    
    errorbar(TF_conc_t,real_data_Ps(:,indv(kk)),real_data_Ps_std(:,indv(kk)),'LineWidth',4)
    hold on
    
    TR = time_dep_TR_new_wSigma(nbd,energyi_v,TF_conc_t,A_conc_t,mut_mat(indv(kk),:),vmax_v);
    plot(TF_conc_t,TR,'+','LineWidth',2)
    xlabel('TF concentration')
    ylabel('transcription rate')
    title(string(title_name(indv(kk))))
end