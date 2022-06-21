function [indi_err,overall_err,overall_TR] = plot_model_data(ii,energyi,TF_conc_t,RNAp_conc_t,vmax_array,group_array,ind)

load('promoter_activity_single.mat')

% load('H_RNAP.mat', 'H_conc_t')
% load('A_RNAP.mat','A_conc_t')

real_data_Ps = Ps_promoter_activity_mean(:,2:9); % exclude the real WT, instead WT is Ps only, no mutations
real_data_Ps_std = Ps_promoter_activity_std(:,2:9);
real_data_Pv = Pv_promoter_activity_mean(:,2:9);
real_data_Pv_std = Pv_promoter_activity_std(:,2:9);


if ii == 'v'
    real_data = real_data_Pv;
    real_data_std = real_data_Pv_std;
    title_name = {'Pv','1*4*','2*4*','3*4*','124*','134*','234*','1234*'};
%     RNAp_conc_t = A_conc_t;
else
    real_data = real_data_Ps;
    real_data_std = real_data_Ps_std;
    title_name = {'Ps','1*','2*','3*','12*','13*','23*','123*'};
%     RNAp_conc_t = H_conc_t;
end

%% assign parameters
nbd = 4;
mut_mat = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[1,0,0];[0,0,0]];
n_strain = length(ind);

%% now see how the newer parameters do
% and record the error for each mutation
indi_err = zeros(length(ind),1);
overall_TR = zeros(length(TF_conc_t),n_strain);
figure();
for kk = 1:length(ind)
    subplot(4,2,kk)
    
    errorbar(TF_conc_t,real_data(:,ind(kk)),real_data_std(:,ind(kk)),'LineStyle','none','LineWidth',2)
    hold on
    
    TR = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNAp_conc_t,mut_mat(ind(kk),:),vmax_array,group_array);
    overall_TR(:,kk) = TR;
    
    indi_err(kk) = weighted_msd(real_data(:,ind(kk)),TR);
    plot(TF_conc_t,TR,'LineWidth',2)
    xlabel('TF concentration')
    ylabel('transcription rate')
    
    ylim([0 1000])
    title(string(title_name(ind(kk))))
end

rd = real_data(:,ind);

overall_err = weighted_msd(rd(:),overall_TR(:));
end