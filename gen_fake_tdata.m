function t_overall = gen_fake_tdata(tspan,nbd,tfunc)
% input: time lasting for the experiment, number of binding site, and how
% concentration varies with time

% output: time dependent configuration profiles ready to be plugged in to
% solve for probability of transcription (note that there will be
% tspan*2^nbd rows of this struct). To calculate PofTR, take every tspan
% number of rows for next step

% total time of experiment 
tall = [0:0.01:tspan];

% number of binding boxes
% nbd = 3;

% list all possible configurations, 
% with last digit being if RNAP bound or not
bins = dec2bin(0:(2^nbd-1), nbd) - '0';

% at each time point, construct the config profile. Assume the number of
% proteins bound stay constant over time. Just want to examine how concentration
% changes gene expression here

t_overall = [];

for tt = 1:tspan % at each time point
    
    for nn = 1:2^nbd % each profile contains 2^n possible configurations, for each configurati
        
        curr_config = bins(nn,:); % print out current configuration
        curr_conc = tfunc(tall(tt));
        [n_arr,c_arr] = gen_prof(curr_config,curr_conc);
        
        % record two fields for each configuration to complete information
        % for one configuration
        cur_config(nn).num_arr = n_arr;
        cur_config(nn).conc_arr = c_arr;
        
    end  
    
    t_overall = [t_overall,cur_config];
    
end



                
    