function [TR_overall_v,TR_overall_s] = time_dep_TR_new(nbd,energyi,TF_conc_t,mut)
%% generate time dependent TR (fake data) to get an intuition about real data
TR_overall_v = zeros(length(TF_conc_t),1);
TR_overall_s = zeros(length(TF_conc_t),1);

for tt = 1:length(TF_conc_t)
    [p_on_v, p_on_s] = transcription_rate_new(nbd,energyi,mut,TF_conc_t(tt));
    TR_overall_v(tt) = p_on_v;
    TR_overall_s(tt) = p_on_s;
end

figure();
plot(TF_conc_t,TR_overall_v,'LineWidth',5)
hold on
plot(TF_conc_t,TR_overall_s,'LineWidth',5)
xlabel('time')
ylabel('transcription rate')
ylim([0 1])
set(gca,'FontSize',17)
legend('Pv','Ps')
end
