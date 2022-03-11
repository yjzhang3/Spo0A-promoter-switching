%% fixed parameters
TF_conc_t = @(x) x^3+30*x;
RNAp_conc = 50;
nbd = 4;
tspan = 2500;
original = [0,20,20,20,0,0,0,0,0,0];
% clone many copies of the original energy array
original_mat = repmat(original,3,1);

%% two TF-RNAp interaction, one activation one repression
close all
change_1 = [400];
pos = 9; % A-RNAp
[TR_all_1,new_E_mat] = change_energy(change_1,pos,RNAp_conc,tspan,TF_conc_t,nbd,original_mat);

change_2 = [-20];
pos = 10; % B-RNAp
[TR_all_2,new_E_mat_2] = change_energy(change_2,pos,RNAp_conc,tspan,TF_conc_t,nbd,new_E_mat);

figure();
for nn = 1:length(TR_all_2(:,1))
    plot(1:tspan,TR_all_2(nn,:),'LineWidth',4)
    ylim([0 1])
    
    spec1 = sprintf('energy config %d',nn);
    yline(TR_all_2(nn,end),'--',spec1)

    xlabel('time')
    ylabel('transcription rate')
    
    hold on
    set(gca,'FontSize',15)
end

