%% a note on raw data
% Ps: column 2 and onward are one promoter Ps with all mutations. First
% column is true wild type (both promoters exist)
% Pv: column 2-9 are one promoter Pv with all mutations. First column is
% true wild type (both promoters intact). Last column is when both promoters are deleted.
% all data start from T2 to T10 (index 1 to 9)
% we provide spo0A dynamics from hour  3 to 8
load('promoter_activity_single.mat')
load('Spo0A_dynamics.mat')
load('Apr21_purePs_estimateRNAp.mat', 'RNAp_conc_t_real')
%% parameters
real_data_Ps = Ps_promoter_activity_mean(2:7,1:8); % exclude the real WT, instead WT is Ps only, no mutations
mut_mat = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[0,0,1];[0,0,0]];

nbd = 5;
TF_conc_t = wtaCopy(1:6:end);
energy_len = nbd+nbd*(nbd-1)/2;

mut_mat = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[0,0,1];[0,0,0]];
lb = zeros(energy_len,1)-20;
ub = zeros(energy_len,1)+15;

lb(1:4) = 5; % standard energy must be positive

% binding box 2 and 4 have lower affinity (high standard energy)
lb(2) = 9;
lb(4) = 6.5; 

lb(5) = 4; % promoter energy should be no greater than standard

 % introduce TF-TF repulsions
lb([6:8 10:11 13]) = 1; 
lb(11) = 5; % box 2 and 4 repulsion is larger
ub([6:8 10:11 13]) = 10; % TF-repulsion shouldn't be dominate

lb(9) = 5; % make one of them necessarily a repressor
ub(9) = 10; % but not a super strong repressor

% normalize the real data
for ll = 1:length(real_data_Ps(1,:))
    real_data_Ps(:,ll) = real_data_Ps(:,ll)./max(real_data_Ps(:,ll));
end

figure();
for kk = 1:length(mut_mat)
    subplot(4,2,kk)
    plot(TF_conc_t,real_data_Ps(:,kk),'o','LineWidth',2)
    ylim([0 1])
    xlabel('TF concentration')
    ylabel('transcription rate')
    title(sprintf('mutation %d',kk))
end

%%
[pars,diff] = fit_data(nbd,TF_conc_t,mut_mat,real_data_Ps,lb,ub);
[M,I] = min(diff);
%% figures

% normalize the real data
for ll = 1:length(real_data_Ps(1,:))
    real_data_Ps(:,ll) = real_data_Ps(:,ll)./max(real_data_Ps(:,ll));
end

figure();
for kk = 1:length(mut_mat)
    subplot(4,2,kk)
    
    plot(TF_conc_t,real_data_Ps(:,kk),'o','LineWidth',4)
    ylim([0 1])
    hold on
    
    TR = time_dep_TR_new(nbd,pars(I,:),TF_conc_t,mut_mat(kk,:));
    plot(TF_conc_t,TR,'LineWidth',4)
    xlabel('TF concentration')
    ylabel('transcription rate')
    title(sprintf('mutation %d',kk))
end