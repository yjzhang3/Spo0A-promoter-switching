function Z = Z_on(nbd,energyi,typei,coopi,TF_conc)
% nbd: number of binding site
% typei: type of each TF (activator > 1, repressor < 1), excluding mRNA
% each number represents the activating or repressive power, which is just
% the factor that accounts for TF-mRNA interaction
% coopi: TF-TF cooperativity of each configuration

% output:
% Z = partition function of all the configurations when transcription
% machinary is bound

%% list all possible configurations, 
% with last digit being if RNAP bound or not
bins = dec2bin(0:(2^nbd-1), nbd) - '0';

% then found configurations that are ON (last digit == 1)
ind = find(bins(:,end)==1);
on_config = bins(ind,:);


%% calculate partition function
Z = 0;

for ll = 1:length(on_config) % for all on configurations
    curr_config = on_config(ll,:); % current configuration, just print out
    
    fe_on = free_energy(curr_config,energyi,TF_conc);
    
    list_int = typei.*curr_config(1:end-1); % calculate TF-mRNA interacton at each site
    % if a binding site is not occupied (0), then it can't have TF-mRNA
    % interaction, which is 0
    
    curr_Q = TF_mRNA_int(list_int,2); % calculate TF-mRNA interaction 
       % influence on the entire configuration (exclude the last digit
       % because it represents mRNA itself)
       % mode == 1, meaning we are using additive model
    
    Z = Z + curr_Q*coopi(ll)*exp(-fe_on);

end
end
            
        
    



