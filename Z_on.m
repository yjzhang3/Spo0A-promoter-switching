function Z = Z_on(nbd,energyi,TF_conc,RNAp_conc,mut)
% nbd: number of binding site
% typei: type of each TF (activator > 1, repressor < 1), excluding mRNA
% each number represents the activating or repressive power, which is just
% the factor that accounts for TF-mRNA interaction
% coopi: TF-TF cooperativity of each configuration
% mut: an array that specifies mutations (should be the same length as the
% number of binding sites - 1, because nbd includes promoter). 

% output:
% Z = partition function of all the configurations when transcription
% machinary is bound

%% list all possible configurations, 
% with last digit being if RNAP bound or not
bins = dec2bin(0:(2^nbd-1), nbd) - '0';

% consider mutated sites
for ss = 1:length(mut)
    if mut(ss) == 0
        bins(:,ss) = 0;
    end
end

% then found configurations that are ON (last digit == 1)
ind = find(bins(:,end)==1);
on_config = bins(ind,:);


%% calculate partition function where mRNA binds
Z = 0;

for ll = 1:length(on_config) % for all on configurations
    curr_config = on_config(ll,:); % current configuration, just print out
    
    Z = Z+Z_perconfig(curr_config,energyi,TF_conc,RNAp_conc);

end

end
            
        
    



