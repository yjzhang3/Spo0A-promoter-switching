%% compare simulated and real data

subplot(1,2,1)
for ll = 1:length(mut_mat)
    TR_overall = time_dep_TR(nbd,x(2:end),TF_conc_t,x(1),mut_mat(:,ll));
    plot(TR_overall,'LineWidth',4)
    xlabel('time (sec)')
    ylabel('transcription rate')
    set(gca,'FontSize',12)
    title('simulated data')
    hold on
end
hold off

subplot(1,2,2)
for ll = 1:length(real_data)
    plot(real_data(:,ll)/max(real_data(:,ll)),'LineWidth',4)
    xlabel('time (hr)')
    ylabel('promoter activity')
    set(gca,'FontSize',12)
    title('experimental data')
    hold on
end
hold off