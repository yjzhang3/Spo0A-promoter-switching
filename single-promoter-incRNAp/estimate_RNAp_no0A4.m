function [RNAp_conc_t,vmax] = estimate_RNAp_no0A4(ii)
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
load('promoter_activity_single.mat')

real_data_Ps = Ps_promoter_activity_mean(:,2:9); % exclude the real WT, instead WT is Ps only, no mutations
real_data_Ps_std = Ps_promoter_activity_std(:,2:9);
real_data_Pv = Pv_promoter_activity_mean(:,2:9);
real_data_Pv_std = Pv_promoter_activity_std(:,2:9);

TF_conc_t = get_aps(2:10);

if ii == 'v'
    real_data = real_data_Pv;
    real_data_std = real_data_Pv_std;
    title_name = {'Pv','1*4*','2*4*','3*4*','124*','134*','234*','1234*'};
    titlen = 'Pv-sigmaA';
end
if ii == 's'
    real_data = real_data_Ps;
    real_data_std = real_data_Ps_std;
    title_name = {'Ps','1*','2*','3*','12*','13*','23*','123*'};
    titlen = 'Ps-sigmaH';
end

%% other parameters

% modify this accrodign to the strain you are fitting to
energyi_sigma = zeros(10,1); % energy arrays used for estimating RNAP concentration
% doesn't matter what these energies are (set promoter energy to be 0, so
% that the found [RNAp] parameter incorporate energy value as well

nbd = 4;
%% fit 123*  data first to find time-dependent RNAP
% treat both energies and RNAp as unknowns

mut_123_data = real_data(:,end); % the 123* strain
mut_123 = [0,0,0];

%% find [RNAp] as functino of time first
[pars,diff] = fit_data_sigma(nbd,energyi_sigma,TF_conc_t,mut_123,mut_123_data);

RNAp_conc_t = pars(1:end-1);
vmax = pars(end);
% the array sigma_conc_t is the time-dependent parameter that incorporates
% [RNAp] and promoter binding affinity. So later when we apply this array
% to all mutation strains, can just let promoter energy be 0 always (keep
% it fixed, because exp(0) = 1 and this will be multiplied by
% [RNAp]-containing terms

%% plot to see the fitting results
figure();
plot(TF_conc_t,mut_123_data,'o','LineWidth',3.5)
hold on
plot(TF_conc_t,time_dep_TR_new_wSigma(nbd,energyi_sigma,TF_conc_t,A_conc_t,mut_123,vmax),'LineWidth',3.5)
xlabel('[Spo0A]')
ylabel('Promoter Activity')
title(titlen)
