function TR_overall = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNAp_conc_t,mut,vmax_array,group_array)
%% generate time dependent TR (fake data) to get an intuition about real data
TR_overall = zeros(length(TF_conc_t),1);

for tt = 1:length(TF_conc_t)
    TR_overall(tt) = transcription_rate_new(nbd,energyi,mut,TF_conc_t(tt),RNAp_conc_t(tt),vmax_array,group_array);
end
% 
% 
% figure();
% plot(TF_conc_t,TR_overall,'LineWidth',5)
% hold on
% xlabel('[Spo0A~P]')
% ylabel('transcription rate')
% % ylim([0 1000])
% set(gca,'FontSize',17)
end
