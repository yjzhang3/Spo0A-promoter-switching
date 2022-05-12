

load('PsPv_promoter_activity.mat')
TF_conc_t = get_aps(2:10);

figure();
for kk = 1:4
    
    errorbar(TF_conc_t,PsPv_avg(:,kk),PsPv_std(:,kk),'LineWidth',4)
    hold on
    title('Wild Type Promoter')
    xlabel('TF Concentration')
    ylabel('Promoter Activity (Miller Unit)')
end
legend('WT','1*','2*','3*')

figure();
for kk = 5:8
    errorbar(TF_conc_t,PsPv_avg(:,kk),PsPv_std(:,kk),'LineWidth',4)
    hold on
    title('Wild Type Promoter')
    xlabel('TF Concentration')
    ylabel('Promoter Activity (Miller Unit)')
end
legend('12*','13*','23*','123*')
