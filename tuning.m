function success = tuning(step,stepsize,pos_all,mRNA_conc,tspan,TF_conc_t,nbd,original)

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
% % A activates mRNA, B repels mRNA
% original = [0.5,0.8,20,0,-5,5];
% % energyi_variation = energyi_normal,power(2*step,num_pos)+1);
% % energyi_variation(:,end) = energyi_normal;


all_combo = nmultichoosek([0-step:0+step],length(pos_all));


for ss = 1:length(all_combo)
    energyi_normal = original;
    for nn = 1:length(pos_all)
        energyi_variation = change_ele(energyi_normal,pos_all(nn),all_combo(ss,nn).*stepsize);
        energyi_normal = energyi_variation;
    end
    
    for tt = 1:tspan
        curr = transcription_rate(nbd,energyi_variation,TF_conc_t(tt),mRNA_conc);
        TR_t_curr(tt) = curr;
    end
    
    plot(1:tspan,TR_t_curr)
    spec1 = sprintf('A-mRNA = %.2f and B-mRNA = %.2f',energyi_variation(pos_all(1)),energyi_variation(pos_all(2))); 
    yline(TR_t_curr(tt-1),'--',spec1)
    xlabel('time')
    ylabel('transcription rate')
    spec2 = sprintf('tuning TF-mRNA attraction (%d sites simultaneously)',length(pos_all));
    title(spec2)
    hold on
    set(gca,'FontSize',15)
end
success = 1;

end
    