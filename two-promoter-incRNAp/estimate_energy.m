function [pars,diff] = estimate_energy(n_strain,ind,file)
%% this file finds the energy for smaller set of data
%% load data and choose strains
% ii tells us Pv or Ps
% n_strain tell us the number of strains to fit with
% ind tells us which strain to fit

load('PsPv_promoter_activity.mat')

load('A_conc_t_reasonable.mat', 'A_conc_t')
load('H_conc_t_reasonable.mat','H_conc_t')

RNApA_conc_t = A_conc_t;
RNApH_conc_t = H_conc_t;

real_data = PsPv_avg;
real_data_std = PsPv_std;

TF_conc_t = get_aps(2:10);

title_name = {'WT','1*','2*','3*','12*','13*','23*','123*'};

titlename = title_name(ind);

%% assign parameters
nbd = 5;
mut_mat = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[0,0,1];[0,0,0]];

%% set up bounds for energy

lb_overall = zeros(39,1)-20;
ub_overall = zeros(39,1)+20;

% lb_overall(1:3) = 1; % binding affinity of each site > 1 (nonzero and greater than promoter affinity)
% lb_overall(7) = 5; % must have one repressor

lb_overall(4:5) = 0;
ub_overall(4:5) = 0; % keep promoter energy fixed because we know it's embedded in concentration terms

% lb_overall(7) = 0;
% ub_overall(7) = 0; % TF bound at 0A1 and 3 probably won't have interaction energy

% vmax should be positive
lb_overall(16:end) = 0;
ub_overall(16:end) = 100000;

% vmax0_ind = 15+[3,4,6,9,12:18,20,21,23,24];
% 
% ub_overall(vmax0_ind) = 0;

%% index smaller dataset
    
mut_mat_new = mut_mat(ind,:); 
real_data_new = real_data(:,ind);

%% start to fit data
[pars,diff] = fit_data(nbd,TF_conc_t,RNApH_conc_t,RNApA_conc_t,mut_mat_new,real_data_new,lb_overall,ub_overall,n_strain,file);

%%
final_energyi = pars(1:15);
final_vmax = pars(16:end);

%% now see how the newer parameters do
figure();
for kk = 1:length(ind)
    subplot(4,2,kk)
%     figure();
    
    errorbar(TF_conc_t,real_data(:,ind(kk)),real_data_std(:,ind(kk)),'LineStyle','none','LineWidth',2)
    hold on

    TR = time_dep_TR_new_wSigma(nbd,final_energyi,TF_conc_t,RNApH_conc_t,RNApA_conc_t,mut_mat(ind(kk),:),final_vmax);
    plot(TF_conc_t,TR,'LineWidth',2)
    xlabel('TF concentration')
    ylabel('transcription rate')
    
    ylim([0 1000])
    title(string(title_name(ind(kk))))
end
% 
% end