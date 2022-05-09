function [pars,diff] = estimate_energy(ii,n_strain,ind)
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

if ii == 'v'
    real_data = real_data_Pv;
    real_data_std = real_data_Pv_std;
    title_name = {'Pv','1*4*','2*4*','3*4*','124*','134*','234*','1234*'};
    RNAp_conc_t = A_conc_t;
else
    real_data = real_data_Ps;
    real_data_std = real_data_Ps_std;
    title_name = {'Ps','1*','2*','3*','12*','13*','23*','123*'};
    RNAp_conc_t = H_conc_t;
end

titlename = title_name(ind);

%% assign parameters
nbd = 4;
mut_mat = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[0,0,1];[0,0,0]];

%% set up bounds for energy

lb_overall = zeros(10,1)-20;
ub_overall = zeros(10,1)+20;

% lb_overall(1:3) = 1; % binding affinity of each site > 1 (nonzero and greater than promoter affinity)
% lb_overall(7) = 5; % must have one repressor

lb_overall(4) = 0;
ub_overall(4) = 0; % keep promoter energy fixed because we know it's embedded in concentration terms

% vmax should be positive
lb_overall(11:end) = 0;
ub_overall(11:end) = 10000;

%% index smaller dataset
    
mut_mat_new = mut_mat(ind,:); 

real_data_new = real_data(:,ind);

%% start to fit data

[pars,diff] = fit_data(nbd,TF_conc_t,RNAp_conc_t,mut_mat_new,real_data_new,lb_overall,ub_overall,n_strain);

%%
final_energyi = pars(1:10);
final_vmax = pars(11:end);

%% now see how the newer parameters do
figure();
for kk = 1:length(ind)
    subplot(4,2,kk)
    
    errorbar(TF_conc_t,real_data(:,ind(kk)),real_data_std(:,ind(kk)),'LineWidth',4)
    hold on
    
    TR = time_dep_TR_new_wSigma(nbd,final_energyi,TF_conc_t,RNAp_conc_t,mut_mat(ind(kk),:),final_vmax);
    plot(TF_conc_t,TR,'+','LineWidth',2)
    xlabel('TF concentration')
    ylabel('transcription rate')
    title(string(title_name(ind(kk))))
end

end
