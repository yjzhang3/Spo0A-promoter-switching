function Z = Z_off(nbd,energyi,TF_conc,RNAp_conc)
% input:
% nbd: number of binding site
% energyi: binidng energy of ea
% coopi: TF-TF cooperativity of each configuration

% output:
% Z = partition function of all the configurations when transcription
% machinary is bound

%% list all possible configurations, 
% with last digit being if RNAP bound or not
bins = dec2bin(0:(2^nbd-1), nbd) - '0';

% then found configurations that are ON (last digit == 1)
ind = find(bins(:,end)==0);
off_config = bins(ind,:);

%% calculate partition function where mRNA does not bind
Z = 0;

for ll = 1:length(off_config) % for all on configurations
    curr_config = off_config(ll,:); % current configuration, just print out
    
    Z = Z+Z_perconfig(curr_config,energyi,TF_conc,RNAp_conc);

end

% figure();
% plot(1:length(off_config),Z_test);
% xlabel('configuration #')
% ylabel('partition function value') 
% title('Z off')

end
            
        
    



