function plot_pars(csvf,perc)

% csvf: the csv file to be read
% perc: up to how many percent away from the minimum error among all sets
% of parameters

M = csvread(csvf);
minerr = min(M(:,end));
uind = (M(:,end)<=minerr*(1+perc));
M_final = M(uind,:);

mu_array = zeros(length(M_final(1,:)),1);
sgm_array = zeros(length(M_final(1,:)),1);
for kk = 1:length(M_final(1,:))
    pd = fitdist(M_final(:,kk),'Normal');
    mu_array(kk) = pd.mu;
    sgm_array(kk) = pd.sigma;
end

%% plot the parameter with their error bar
% first plot the energy (we know it must be the first 10)
figure();
bar(mu_array(1:10));
hold on
errorbar(mu_array(1:10),sgm_array(1:10),'LineStyle','none','LineWidth',2);
xlabel('energy types')
xticklabels({'0A1','0A2','0A3','p','12','13','1P','23','2p','3p'})
hold off
saveas(gcf,append(csvf(1:end-4),'-energy','.jpeg'))

%% plot the parameter with their error bar
% then plot the maximum initiation rate
figure();
bar(mu_array(11:end-1));
hold on
errorbar(mu_array(11:end-1),sgm_array(11:end-1),'LineStyle','none','LineWidth',2);
xlabel('v_max')
saveas(gcf,append(csvf(1:end-4),'-vmax','.jpeg'))
