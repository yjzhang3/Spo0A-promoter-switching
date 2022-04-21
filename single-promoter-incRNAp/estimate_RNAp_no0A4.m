%% this file is trying to find how RNAp (more specificlaly, sigma-H bound
% holonenzyme concentration increaeses with time 
load('promoter_activity_single.mat')
load('Spo0A_dynamics.mat')
real_data_Ps = Ps_promoter_activity_mean(2:7,1:8); % exclude the real WT, instead WT is Ps only, no mutations

TF_conc_t = wtaCopy(1:6:end);
RNAp_conc_t_mock = zeros(length(TF_conc_t),1)+1;

nbd = 4;
mut_mat = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[0,0,1];[0,0,0]];

%% re-watch what is the issue of previous fitting
% use the energy parameters found for the previous best run
load('Apr13_results_purePs_onerep_3TFTFrepsul.mat', 'pars');
energyi = pars(1,:);

for ll = 1:length(real_data_Ps(1,:))
    real_data_Ps(:,ll) = real_data_Ps(:,ll)./max(real_data_Ps(:,ll));
end

figure();
for kk = 1:length(real_data_Ps)
    subplot(4,2,kk)
    
    plot(TF_conc_t,real_data_Ps(:,kk),'o','LineWidth',4)
    ylim([0 1])
    hold on
    
    TR = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNAp_conc_t_mock,mut_mat(kk,:));
    plot(TF_conc_t,TR,'LineWidth',4)
    xlabel('TF concentration')
    ylabel('transcription rate')
    title(sprintf('mutation %d',kk))
end
%% fit 123*  data
% treat both eenergies and RNAp as unknowns
load('Apr13_results_purePs_onerep_3TFTFrepsul.mat', 'lb');
load('Apr13_results_purePs_onerep_3TFTFrepsul.mat', 'ub');

lb_overall = zeros(16,1);
ub_overall = zeros(16,1);
lb_overall(1:10) = lb;
ub_overall(1:10) = ub;

lb_overall(11:end) = zeros(length(TF_conc_t),1);
ub_overall(11:end) = zeros(length(TF_conc_t),1)+1000000;

all_mut_data = real_data_Ps(:,end); % the 123* strain
mut_RNAp = [0,0,0];

%%
[new_pars,diff] = fit_data_sigma(nbd,TF_conc_t,mut_RNAp,all_mut_data,lb_overall,ub_overall);

%% now see  how the newe parameters do
RNAp_conc_t_real = new_pars(11:end);
energyi = new_pars(1:10);

figure();
for kk = 1:length(real_data_Ps)
    subplot(4,2,kk)
    
    plot(TF_conc_t,real_data_Ps(:,kk),'o','LineWidth',4)
    ylim([0 1])
    hold on
    
    TR = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNAp_conc_t_real,mut_mat(kk,:));
    plot(TF_conc_t,TR,'LineWidth',4)
    xlabel('TF concentration')
    ylabel('transcription rate')
    title(sprintf('mutation %d',kk))
end