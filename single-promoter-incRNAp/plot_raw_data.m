ii = 's'
real_data = [];
real_data_std = [];

load('promoter_activity_single.mat')
load('Spo0A_dynamics.mat')

real_data_Ps = Ps_promoter_activity_mean(2:7,1:8); % exclude the real WT, instead WT is Ps only, no mutations
real_data_Ps_std = Ps_promoter_activity_std(2:7,1:8);
real_data_Pv = Pv_promoter_activity_mean(2:7,1:8);
real_data_Pv_std = Ps_promoter_activity_std(2:7,1:8);

wta = wtaCopy(1:6:end);

if ii == 'v'
    real_data = real_data_Pv;
    real_data_std = real_data_Pv_std;
else
    real_data = real_data_Ps;
    real_data_std = real_data_Ps_std;
end

% title_name = {'Ps','1*','2*','3*','12*','13*','23*','123*'};
figure();
for kk = 5:length(real_data)
%     subplot(4,1,kk)
    
    errorbar(wta,real_data(:,kk),real_data_std(:,kk),'LineWidth',4)
    hold on
    
    xlabel('TF Concentration')
    ylabel('Promoter Activity (Miller Unit)')
%     title(string(title_name(kk)))
end
legend('12*','13*','23*','123*')
