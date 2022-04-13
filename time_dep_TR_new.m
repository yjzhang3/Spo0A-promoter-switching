function TR_overall = time_dep_TR_new(nbd,energyi,TF_conc_t,mut)
%% generate time dependent TR (fake data) to get an intuition about real data
TR_overall = zeros(length(TF_conc_t),1);

for tt = 1:length(TF_conc_t)
    TR_overall(tt) = transcription_rate_new(nbd,energyi,mut,TF_conc_t(tt));
end

% figure();
% plot(1:tspan,TR_overall,'LineWidth',5)
% hold on
% xlabel('time')
% ylabel('transcription rate')
% ylim([0 1])
% set(gca,'FontSize',17)
end
