function plot_raw_data(ii)

load('promoter_activity_single.mat')
load('Spo0A_dynamics.mat')

real_data_Ps = Ps_promoter_activity_mean(2:7,2:9); % exclude the real WT, instead WT is Ps only, no mutations
real_data_Ps_std = Ps_promoter_activity_std(2:7,2:9);
real_data_Pv = Pv_promoter_activity_mean(2:7,2:9);
real_data_Pv_std = Pv_promoter_activity_std(2:7,2:9);

wta = wtaCopy(1:6:end);

if ii == 'v'
    real_data = real_data_Pv;
    real_data_std = real_data_Pv_std;
    titlen = 'Pv-sigmaA';
end
if ii == 's'
    real_data = real_data_Ps;
    real_data_std = real_data_Ps_std;
    titlen = 'Ps-sigmaH';
end

figure();
for kk = 1:4
    
    errorbar(wta,real_data(:,kk),real_data_std(:,kk),'LineWidth',4)
    hold on
    title(titlen)
    xlabel('TF Concentration')
    ylabel('Promoter Activity (Miller Unit)')
end
legend('WT','1*','2*','3*')

figure();
for kk = 5:length(real_data)
    errorbar(wta,real_data(:,kk),real_data_std(:,kk),'LineWidth',4)
    hold on
    title(titlen)
    xlabel('TF Concentration')
    ylabel('Promoter Activity (Miller Unit)')
end
legend('12*','13*','23*','123*')
