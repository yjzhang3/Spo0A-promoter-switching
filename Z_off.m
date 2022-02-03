function Z = Z_off(nbd,config_prof,G0)
% input:
% nbd: number of binding site
% config_prof: used for calculating the free energy of each configuration.
% It'll be a struct array, and each element has two fields, conc_arr and
% num_arr

% output:
% Z = partition function of all the configurations when transcription
% machinary is bound

%% list all possible configurations, 
% with last digit being if RNAP bound or not
bins = dec2bin(0:(2^nbd-1), nbd) - '0';

% then found configurations that are ON (last digit == 1)
ind = find(bins(:,end)==0);
off_config = bins(ind,:);

%% calculate partition function
Z = 0;
Z_test = zeros(length(off_config),1);

for ll = 1:length(off_config) % for all on configurations
    curr_config = off_config(ll,:); % current configuration, just print out
    
    % index the conc and num profile of this configuration
    conc_arr = config_prof(ind(ll)).conc_arr;
    num_arr = config_prof(ind(ll)).num_arr;
    
    fe_off = free_energy(G0,conc_arr,num_arr);
    Z = Z + exp(-fe_off);
    Z_test(ll) = Z;
end

% figure();
% plot(1:length(off_config),Z_test);
% xlabel('configuration #')
% ylabel('partition function value') 
% title('Z off')

end
            
        
    



