%% this file is trying to find how RNAp (more specificlaly, sigma-H bound
% holonenzyme concentration increaeses with time), and then incorporate
% this inofmration back to other strain types

% Ps: column 2 and onward are one promoter Ps with all mutations. First
% column is true wild type (both promoters exist)
% Pv: column 2-9 are one promoter Pv with all mutations. First column is
% true wild type (both promoters intact). Last column is when both promoters are deleted.
% all data start from T2 to T10 (index 1 to 9)
% we provide spo0A dynamics from hour  3 to 8 (so index into 2 to 7)

% no 0A4 box for now
ii = 'v'
real_data = [];
real_data_std = [];

load('promoter_activity_single.mat')
load('Spo0A_dynamics.mat')

real_data_Ps = Ps_promoter_activity_mean(2:7,1:8); % exclude the real WT, instead WT is Ps only, no mutations
real_data_Ps_std = Ps_promoter_activity_std(2:7,1:8);
real_data_Pv = Pv_promoter_activity_mean(2:7,1:8);
real_data_Pv_std = Ps_promoter_activity_std(2:7,1:8);

wta = wtaCopy(1:6:end);

if ii == 'v'
    real_data = real_data_Pv;
    real_data_std = real_data_Pv_std;
else
    real_data = real_data_Ps;
    real_data_std = real_data_Ps_std;
end
%% other parameters

% modify this accrodign to the strain you are fitting to
TF_conc_t = wta;
energyi_sigma = zeros(10,1); % energy arrays used for estimating RNAP concentration
% doesn't matter what these energies are (set promoter energy to be 1, so
% that the found [RNAp] parameter incorporate energy value as well

nbd = 4;
mut_mat = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[0,0,1];[0,0,0]];

%% fit 123*  data first to find time-dependent RNAP
% treat both energies and RNAp as unknowns

mut_123_data = real_data(:,end); % the 123* strain
mut_123 = [0,0,0];

%% find [RNAp] as functino of time first
[pars_1,diff] = fit_data_sigma(nbd,energyi_sigma,TF_conc_t,mut_123,mut_123_data);

sigmaH_conc_t = pars_1(1:end-1);
vmax_123 = pars_1(end);
% the array sigma_conc_t is the time-dependent parameter that incorporates
% [RNAp] and promoter binding affinity. So later when we apply this array
% to all mutation strains, can just let promoter energy be 0 always (keep
% it fixed, because exp(0) = 1 and this will be multiplied by
% [RNAp]-containing terms

%% plot to see the fitting results
figure();
plot(TF_conc_t,mut_123_data,'o','LineWidth',3.5)
hold on
plot(TF_conc_t,time_dep_TR_new_wSigma(nbd,energyi_sigma,TF_conc_t,sigmaH_conc_t,mut_123,vmax_123),'LineWidth',3.5)
xlabel('[Spo0A]')
ylabel('Promoter Activity')

%% then find energy (bounds setup)

lb_overall = zeros(11,1)-20;
ub_overall = zeros(11,1)+20;

% lb_overall(1:3) = 1; % binding affinity of each site > 1 (nonzero and greater than promoter affinity)
% lb_overall(7) = 5; % must have one repressor

lb_overall(4) = 0;
ub_overall(4) = 0; % keep promoter energy fixed because we know it's embedded in concentration terms

% vmax should be positive
lb_overall(end) = 0;
ub_overall(end) = 10000;

%% fit all 8 strains together
[pars_2,diff] = fit_data_new(nbd,TF_conc_t,sigmaH_conc_t,mut_mat,real_data,lb_overall,ub_overall);
final_energyi = pars_2(1:10);
final_vmax = pars_2(11);

%% now see how the newer parameters do
title_name = {'Ps','1*','2*','3*','12*','13*','23*','123*'};
figure();
for kk = 1:length(real_data)
    subplot(4,2,kk)
    
    errorbar(TF_conc_t,real_data(:,kk),real_data_std(:,kk),'LineWidth',4)
    hold on
    
    TR = time_dep_TR_new_wSigma(nbd,final_energyi,TF_conc_t,sigmaH_conc_t,mut_mat(kk,:),final_vmax);
    plot(TF_conc_t,TR,'+','LineWidth',2)
    xlabel('TF concentration')
    ylabel('transcription rate')
    title(string(title_name(kk)))
end