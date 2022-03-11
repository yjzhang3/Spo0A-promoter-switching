function TR_overall = time_dep_TR(TF_conc_t,tspan,RNAp_conc,nbd,energyi)
%% generate time dependent TR (fake data) to get an intuition about real data
TR_overall = zeros(tspan,1);

for tt = 0:tspan-1
    TR_overall(tt+1) = transcription_rate(nbd,energyi,TF_conc_t(tt),RNAp_conc);
end

% figure();
% plot(1:tspan,TR_overall,'LineWidth',5)
% hold on
% xlabel('time')
% ylabel('transcription rate')
% ylim([0 1])
% set(gca,'FontSize',17)
end
