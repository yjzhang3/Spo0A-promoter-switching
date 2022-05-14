load('A_conc_t_reasonable.mat', 'A_conc_t')
load('H_conc_t_reasonable.mat', 'H_conc_t')

load('PsPv_promoter_activity.mat')
WT_data = PsPv_avg(:,1);

%%

nbd = 5;
% tentitatively using the energy estimated from May14 (I forgot to set
% promoter energy of Pv to be 0!!)

pars = [20.000000,-3.634995,7.309286,0.000000,5.460547,-14.012441,0.000000,19.999953,-7.793040,-19.980537,8.050116,-4.554194,-7.816389,0.163034,20.000000,62406.969859,823.376199,99998.864091,84812.676014,860.564211,0.013123,947.878427,67605.014190,99787.776563,544.108139,510.822735,0.152240,99914.162715,796.266275,99773.838918,0.080250,98786.453993,3401.428143,551.578026,99842.550686,98666.895079,458.705978,96875.172503,0.041337];

energyi = pars(1:15);
vmax = pars(16:end);

TF_conc_t = get_aps(2:10);
mut_mat = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[0,0,1];[0,0,0]];

%%
time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,H_conc_t,A_conc_t,mut_mat(1,:),vmax);
% hold on 
% scatter(TF_conc_t,WT_data)

%%
mut = mut_mat(1,:);
on_config = dec2bin(0:(2^nbd-1), nbd) - '0';
test = sum(on_config(:,4:5)');
on_ind = find(test > 0);
on_config = on_config(on_ind,:); 

%% plot probabilities of each configuration

figure();
for cc = 1:length(on_config)
    config = on_config(cc,:);
    
    prob_t = zeros(length(TF_conc_t),1);
    for tt = 1:length(TF_conc_t)
        prob_t(tt) = prob_per_config_new(nbd, config,energyi,mut,TF_conc_t(tt),H_conc_t(tt),A_conc_t(tt));
    end
    
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

figure();
for cc = 1:length(on_config)
    config = on_config(cc,:);
    
    prob_t = zeros(length(TF_conc_t),1);
    for tt = 1:length(TF_conc_t)
        prob_t(tt) = prob_per_config_new(nbd, config,energyi,mut,TF_conc_t(tt),H_conc_t(tt),A_conc_t(tt));
    end
    TR_t = prob_t*vmax(cc);
    
    subplot(6,4,cc)
    plot(TF_conc_t,TR_t,'LineWidth',3.5)
    xlabel('[Spo0A]')
    ylabel('Transcription rate')
    yline(max(TR_t),'-',sprintf('p = %.2f',max(TR_t)));
    ylim([0 500])
    
    str = string(config);
    curr_config_name = append(str(1),str(2),str(3),str(4),str(5));
    title(convertStringsToChars(curr_config_name))
end

%% extract those with non-zero probability
key = [[0,1,1,0,1];[0,1,1,1,0];[1,1,1,0,1]];
vmax_key = vmax([10,11,22]);

for cc = 1:3
    config = key(cc,:);
    
    prob_t = zeros(length(TF_conc_t),1);
    for tt = 1:length(TF_conc_t)
        prob_t(tt) = prob_per_config_new(nbd, config,energyi,mut,TF_conc_t(tt),H_conc_t(tt),A_conc_t(tt));
    end
    TR_t = prob_t*vmax_key(cc);
    
    plot(TF_conc_t,TR_t,'LineWidth',3.5)
    hold on
    
    xlabel('[Spo0A]')
    ylabel('Transcription rate')
end
legend('01101','01110','11101')
        
        
        
        