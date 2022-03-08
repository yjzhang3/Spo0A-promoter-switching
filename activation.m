%% fixed parameters
TF_conc_t = @(x) x^3+30*x;
RNAp_conc = 50000;
nbd = 4;
tspan = 400;
original = [10,10,10,10,0,0,0,0,0,0];
% clone many copies of the original energy array
original_mat = repmat(original,5,1);

%% repression from single TF-RNAp interaction
change = [0,-2,-4,-8,-10];
pos = 7; % A-RNAp
[TR_all,~] = change_energy(change,pos,RNAp_conc,tspan,TF_conc_t,nbd,original_mat);

for nn = 1:length(TR_all(:,1))
    plot(1:tspan,TR_all(nn,:),'LineWidth',4)
    
    spec1 = sprintf('energy config %d',nn);
    yline(TR_all(nn,end),'--',spec1)

    xlabel('time')
    ylabel('transcription rate')
    
    hold on
    set(gca,'FontSize',15)
end


%% repression from two TF-RNAp interaction
close all
change = [0,-2,-4,-8,-10];
pos = 7; % A-RNAp
[~,new_E_mat] = change_energy(change,pos,RNAp_conc,tspan,TF_conc_t,nbd,original_mat);

change = [0,-2,-4,-8,-10];
pos = 9; % B-RNAp
[TR_all_2,~] = change_energy(change,pos,RNAp_conc,tspan,TF_conc_t,nbd,new_E_mat);

figure();
for nn = 1:length(TR_all_2(:,1))
    plot(1:tspan,TR_all_2(nn,:),'LineWidth',4)
    
    spec1 = sprintf('energy config %d',nn);
    yline(TR_all_2(nn,end),'--',spec1)

    xlabel('time')
    ylabel('transcription rate')
    
    hold on
    set(gca,'FontSize',15)
end

