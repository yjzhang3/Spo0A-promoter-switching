function TR_overall = time_dep_TR(TF_conc_t,tspan,mRNA_conc,nbd,energyi_multi)
%% generate time dependent TR (fake data) to get an intuition about real data
TR_overall = zeros(tspan,7);

for mm = 1:length(energyi_multi(:,1))
    for tt = 1:tspan
        TR_overall(mm,tt) = transcription_rate(nbd,energyi_multi(mm,:),TF_conc_t(tt),mRNA_conc);
    end
    plot(1:tspan,TR_overall(mm,:))
    hold on
    xlabel('time')
    ylabel('transcription rate')
    set(gca,'FontSize',17)
end
