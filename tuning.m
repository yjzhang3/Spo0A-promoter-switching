function [E_all,TR_all] = tuning(step,stepsize,pos_all,RNAp_conc,tspan,TF_conc_t,nbd,original)

% input:
% step = 2; % how many steps we want to test up and down
% % then we would generate 2*step more energy arrays
% % total there is 2*step+1 arrays
% stepsize = 1.5; % how much energy per step are you tuning?
% 
% pos_all = [5,6]; % which energy value are you tuning?
% 
% TF_conc_t = @(t) 2*t;
% tspan = 10;
% nbd = 3;
% mRNA_conc = 1000;
% 
% original: original energy array

% output: new transcription rate


% generate all possible energy values at the position we want to tune 
all_combo = nmultichoosek([0-step:0+step],length(pos_all));


for ss = 1:length(all_combo) % for each variation 
    energyi_normal = original;
    
    % generate a different energy array
    for nn = 1:length(pos_all) % change each position
        energyi_variation = change_ele(energyi_normal,pos_all(nn),all_combo(ss,nn).*stepsize);
        energyi_normal = energyi_variation;
    end
    
    for tt = 1:tspan % then plot time-dependent TR
        curr = transcription_rate(nbd,energyi_variation,TF_conc_t(tt),RNAp_conc);
        TR_t_curr(tt) = curr;
    end
    
    E_all(ss,:) = energyi_variation;
    TR_all(ss,:) = TR_t_curr;
    
    plot(1:tspan,TR_t_curr,'LineWidth',4)
    spec1 = sprintf('energy config %d',ss);
    yline(TR_t_curr(tt-1),'--',spec1)
    xlabel('time')
    ylabel('transcription rate')
    spec2 = sprintf('tuning p-p attraction (%d energy values simultaneously)',length(pos_all));
    title(spec2)
    hold on
    set(gca,'FontSize',15)
    
end


end
    