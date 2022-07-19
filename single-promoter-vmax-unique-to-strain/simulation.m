%% overall transcription rate
clear all
nbd = 4;

load('promoter_activity_single.mat')

load('H_RNAP.mat', 'H_conc_t')
load('A_RNAP.mat','A_conc_t')

real_data_Ps = Ps_promoter_activity_mean(:,2:9); % exclude the real WT, instead WT is Ps only, no mutations
real_data_Ps_std = Ps_promoter_activity_std(:,2:9);
real_data_Pv = Pv_promoter_activity_mean(:,2:9);
real_data_Pv_std = Pv_promoter_activity_std(:,2:9);

% % Pv best
% load('Jun30-Pv-2group.mat')
% load('A_RNAP.mat','A_conc_t')
% mut = [1,1,1];
% RNAp_conc_t = A_conc_t;
% TF_conc_t = get_aps(2:10);

% Ps best
load('Jun30-Ps-4group.mat')
load('H_RNAP.mat','H_conc_t')
mut = [1,1,1];
RNAp_conc_t = A_conc_t;
TF_conc_t = get_aps(2:10);

energyi = final_energyi;
TR_overall = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNAp_conc_t,mut,vmax_raw(1));
figure();
plot(TF_conc_t,TR_overall,'LineWidth',4)
hold on
errorbar(TF_conc_t,real_data_Pv(:,1),real_data_Pv_std(:,1),'LineStyle','none','LineWidth',2)
xlabel('[Spo0A~P]')
ylabel('overall TR')
ylim([0 300])
set(gca,'FontSize',12)
%% individual probabilities 

bins = dec2bin(0:(2^nbd-1), nbd) - '0';
on_ind = find(bins(:,end) == 1); 
on_config = bins(on_ind,:);
final_config = on_config(ind,:);

figure();
for ii = 1:length(ind)
    curr_ind = ind(ii); % now get the actual index with respect to 8 configurations
    curr_config = final_config(ii,:);
    subplot(4,2,ii)
    
    curr_p = zeros(length(TF_conc_t),1);
    
    for kk = 1:length(TF_conc_t)
        curr_p(kk) = prob_per_config_new(nbd, curr_config, energyi, mut,TF_conc_t(kk),RNAp_conc_t(kk));
    end
    
    curr_vmax = vmax_raw(1); % for strain depdnent vmax, every configuraiton in wild type
    % should have the same vmax (group1)
    
    curr_TR = curr_vmax*curr_p;
    
    
    plot(TF_conc_t,curr_TR,'LineWidth',3)
    
    xlabel('[TF]')
    ylabel('TR')
    ylim([0 200])
    title(append(dec2bin(curr_ind-1,3),'1'))
    set(gca,'FontSize',12)
end