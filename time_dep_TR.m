function TR_overall = time_dep_TR(TF_conc_t,tspan,RNAp,nbd,energyi)
%% generate time dependent TR (fake data) to get an intuition about real data
TR_overall = zeros(tspan,1);

for tt = 1:tspan
    TR_overall(tt) = transcription_rate(nbd,energyi,TF_conc_t(tt),RNAp);
end

% figure();
plot(1:tspan,TR_overall,'LineWidth',5)
xlabel('time')
ylabel('transcription rate')
set(gca,'FontSize',17)
end
