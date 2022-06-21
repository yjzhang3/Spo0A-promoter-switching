%% Pv only

% batch 1: fit to all, but separate 000 and 100 by using different vmax for
% them
clear all
n_strain = 8;
ind = 1:8;
ii = 'v';
vmax_array = struct();
group_array.g1 = [1:4,6:8];
group_array.g2 = 5;
file = 'Jun8-Pvonly-2vmax-b1.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('Jun8-Pvonly-2vmax-b1.mat')

%% Pv only

% batch 2: fit to 23* only and see if separating 000 and 100 by using different vmax for
% them works
clear all
n_strain = 1;
ind = 7;
ii = 'v';
vmax_array = struct();
group_array.g1 = 1;
group_array.g2 = [2:8];
file = 'Jun8-Pvonly-2vmax-b2.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('Jun8-Pvonly-2vmax-b2.mat')

% now the fit is better

%% batch 3: can we improve better? this is only fit one set of data! should be good
clear all
n_strain = 1;
ind = 7;
ii = 'v';
vmax_array = struct();
group_array.g1 = 1;
group_array.g2 = 5;
group_array.g3 = [2:4,6:8];
file = 'Jun8-Pvonly-3vmax-b3.txt';

[bounds,final_energyi,vmax_array,diff] = estimate_energy(ii,n_strain,ind,file,group_array,vmax_array);
save('Jun8-Pvonly-3vmax-b3.mat')

% now the fit is only slightly better. Why can't we fit to 234* only well??
% I tried two groups, 1 and the rest or 4 and the rest, or 1, 4, and the
% rest, but none of them work. What if I try one universal vmax?


%% let's plot probaiblity and TR of 000 and 111 in strain 23*

clear all
close all
nbd = 4;
load('Jun8-Pvonly-3vmax-b3.mat')
load('A_RNAP.mat')
energyi = final_energyi;

% only energy values matter are 1,p,1p
TF_conc_t = get_aps(2:10);
RNAp_conc_t = A_conc_t;
config1 = [0,0,0,1];
config2 = [1,0,0,1];

energyi([2,3,4,5,6,8,9,10]) = 0;
energyi(7) = -4;

RNAp_conc_t(6) = 15;
% TF_conc_t(6) = 0.8;

vmax_array.g1 = 200;


mut = [1,0,0];

p1 = zeros(length(TF_conc_t),1);
p2 = zeros(length(TF_conc_t),1);

for kk = 1:length(TF_conc_t)
    p1(kk) = prob_per_config_new(nbd, config1, energyi,mut,TF_conc_t(kk),RNAp_conc_t(kk));
    p2(kk) = prob_per_config_new(nbd, config2, energyi,mut,TF_conc_t(kk),RNAp_conc_t(kk));
end

TR1 = p1*vmax_array.g1;
TR2 = p2*vmax_array.g2;

figure();
plot(TF_conc_t,p1,'LineWidth',3)
hold on
plot(TF_conc_t,p2,'LineWidth',3)
ylim([0 1])
xlabel('TF concentration')
ylabel('probability')
legend('000','100')

figure();
plot(TF_conc_t,TR1,'LineWidth',3)
hold on
plot(TF_conc_t,TR2,'LineWidth',3)
xlabel('TF concentration')
ylabel('TR')
legend('000','100')

for tt = 1:length(TF_conc_t)
    TR_overall(tt) = transcription_rate_new(nbd,energyi,mut,TF_conc_t(tt),RNAp_conc_t(tt),vmax_array,group_array);
end
% 

% plot the original 23*
[indi_err,overall_err,overall_TR] = plot_model_data(ii,energyi,TF_conc_t,RNAp_conc_t,vmax_array,group_array,ind);

