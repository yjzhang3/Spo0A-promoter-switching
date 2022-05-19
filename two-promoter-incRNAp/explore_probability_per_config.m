load('A_conc_t_reasonable.mat', 'A_conc_t')
load('H_conc_t_reasonable.mat', 'H_conc_t')

load('PsPv_promoter_activity.mat')
WT_data = PsPv_avg(:,1);
WT_data_std = PsPv_std(:,1);

%%
nbd = 5;
% tentitatively using the energy estimated from May14 (I forgot to set
% promoter energy of Pv to be 0!!)

% % this is from first txt file
% pars = [20.000000,-3.634995,7.309286,0.000000,5.460547,-14.012441,0.000000,19.999953,-7.793040,-19.980537,8.050116,-4.554194,-7.816389,0.163034,20.000000,62406.969859,823.376199,99998.864091,84812.676014,860.564211,0.013123,947.878427,67605.014190,99787.776563,544.108139,510.822735,0.152240,99914.162715,796.266275,99773.838918,0.080250,98786.453993,3401.428143,551.578026,99842.550686,98666.895079,458.705978,96875.172503,0.041337];
% 
% % this is from 3rd txt file
% % pars = [19.917175,19.998902,20.000000,0.000000,0.000000,3.805853,0.000000,-17.602240,-6.745312,-19.999993,-6.878131,-16.802702,-19.997049,15.399456,10.006812,384.204469,1124.951660,0.000000,0.000000,1085.502807,0.000000,3249.178580,254.343756,0.000000,616.626749,453.639267,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,99823.026607,0.000000,0.000000,79322.394867,0.000000,0.000000];

% this is from no constraint 
% pars = [8.071011,11.633282,19.976147,0.000000,0.000000,-20.000000,19.103909,-12.932769,1.665981,-19.999906,19.999444,-0.160786,6.309552,-18.297504,-1.176370,331.275547,0.000000,915.695980,1734.513786,93879.104560,90130.600201,99842.471864,22929.408696,7262.405242,433.715751,1978.382269,99192.687733,649.358790,0.000000,758.874026,12536.598074,1.880518,86615.659992,13772.180582,17.398641,99998.879461,46105.913152,19950.182552,1363.363599];
pars = [8.071011,11.633282,19.976147,0.000000,0.000000,-20.000000,19.103909,-12.932769,1.665981,-19.999906,19.999444,-0.160786,6.309552,-18.297504,-1.176370,331.275547,0.000000,915.695980,1734.513786,93879.104560,90130.600201,99842.471864,22929.408696,7262.405242,433.715751,1978.382269,99192.687733,649.358790,0.000000,758.874026,12536.598074,1.880518,86615.659992,13772.180582,17.398641,99998.879461,46105.913152,19950.182552,1363.363599];

% this is from May17 no constraint

energyi = pars(1:15);
vmax = pars(16:end);


TF_conc_t = get_aps(2:10);
mut_mat = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[0,0,1];[0,0,0]];

%% this is only to test if, given a set of optimized energies, can vamx be zero for
% some configurations but still lead to the same TR
% vmax0_ind = [3,4,6,9,12:18,20,21,23,24];
% vmax(vmax0_ind) = 0;
%% compare overall TR between real and simulated
figure();
time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,H_conc_t,A_conc_t,mut_mat(1,:),vmax);
hold on 
errorbar(TF_conc_t,WT_data,WT_data_std,'LineStyle','none','LineWidth',2)

%% list all possible configurations
mut = mut_mat(1,:);
on_config = dec2bin(0:(2^nbd-1), nbd) - '0';
test = sum(on_config(:,4:5)');
on_ind = find(test > 0);
on_config = on_config(on_ind,:); 

%% plot probabilities of each configuration
% 24 configurations (24 rows), 9 time points (9 columns)

figure();
for cc = 1:length(on_config)
    config = on_config(cc,:);
    
    prob_t = zeros(length(TF_conc_t),1);
    for tt = 1:length(TF_conc_t)
        prob_t(tt) = prob_per_config_new(nbd, config,energyi,mut,TF_conc_t(tt),H_conc_t(tt),A_conc_t(tt));
    end
    prob_t_all(cc,:) = prob_t;
    
    subplot(6,4,cc)
    plot(TF_conc_t,prob_t,'LineWidth',3.5)
    xlabel('[Spo0A]')
    ylabel('Probability')
    ylim([0 1])
    yline(max(prob_t),'-',sprintf('p = %.2f',max(prob_t)));
    
    str = string(config);
    curr_config_name = append(str(1),str(2),str(3),str(4),str(5));
    title(convertStringsToChars(curr_config_name))
end
    
%% plot transcription rate of each configuration (prob * vmax)   
% 24 configurations (24 rows), 9 time points (9 columns)

figure();
for cc = 1:length(on_config)
    config = on_config(cc,:);
    
    prob_t = zeros(length(TF_conc_t),1);
    for tt = 1:length(TF_conc_t)
        prob_t(tt) = prob_per_config_new(nbd, config,energyi,mut,TF_conc_t(tt),H_conc_t(tt),A_conc_t(tt));
    end
    TR_t = prob_t*vmax(cc);
    TR_t_all(cc,:) = TR_t;
    
    subplot(6,4,cc)
    plot(TF_conc_t,TR_t,'LineWidth',3.5)
    xlabel('[Spo0A]')
    ylabel('Transcription rate')
    yline(max(TR_t),'-',sprintf('TR = %.4f',max(TR_t)));
    ylim([0 1000])
    
    str = string(config);
    curr_config_name = append(str(1),str(2),str(3),str(4),str(5));
    title(convertStringsToChars(curr_config_name))
end

%% extract those with non-zero probability
mean_prob = mean(prob_t_all');
mean_TR  = mean(TR_t_all');
ind0_prob = find(mean_prob > 0.1);
ind0_TR = find(mean_TR > 100);
ind0_final = intersect(ind0_prob,ind0_TR);

%%
% ind0_final = [3,10,15];
key = on_config(ind0_final,:);
vmax_key = vmax(ind0_final);

%% generate legend
% 
% for gg = 1:length(vmax_key)
%     str = string(on_config(gg,:));
%     curr_config = append(str(1),str(2),str(3),str(4),str(5));
%     legend_list_final(gg) = curr_config;
% end

%%
figure();
for cc = 1:length(vmax_key)
    config = key(cc,:);
%     subplot(2,4,cc)
    
    str = string(config);
    str_im = append(str(1),str(2),str(3),str(4),str(5));
    name = convertStringsToChars(str_im);
    
    prob_t = zeros(length(TF_conc_t),1);
    for tt = 1:length(TF_conc_t)
        prob_t(tt) = prob_per_config_new(nbd, config,energyi,mut,TF_conc_t(tt),H_conc_t(tt),A_conc_t(tt));
    end
    TR_t = prob_t*vmax_key(cc);
    
    plot(TF_conc_t,TR_t,'LineWidth',3.5,'DisplayName',name)
    ylim([0 800])
    
    legend('-DynamicLegend');
    legend('show');
    drawnow;
    
    hold on
    
    xlabel('[Spo0A]')
    ylabel('Transcription rate')
    title('Configurations with top transcription rates')
end
legend(legend_list)        

    
        