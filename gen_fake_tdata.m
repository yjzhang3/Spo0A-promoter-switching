function t_overall = gen_fake_tdata(tspan,nbd,tfunc)
% input: time lasting for the experiment, number of binding site, and how
% concentration varies with time

% output: time dependent configuration profiles ready to be plugged in to
% solve for probability of transcription (note that there will be
% tspan*2^nbd rows of this struct). To calculate PofTR, take every tspan
% number of rows for next step

% total time of experiment 
% say that TF concentration changes with time linearly 
tspan = 10;

% number of binding boxes
nbd = 3;

% list all possible configurations, 
% with last digit being if RNAP bound or not
bins = dec2bin(0:(2^nbd-1), nbd) - '0';

% initialize the array representing the number of proteins bound at each
% binding box. num_arr(i) = j means at ith binding box, there are j
% proteins bound. length of num_arr = number of binding boxes
num_arr = zeros(nbd,1); 
conc_arr = zeros(nbd,1); 

% at each time point, construct the config profile. Assume the number of
% proteins bound stay constant over time. Just want to examine how concentration
% changes gene expression here

t_overall = [];

for tt = 1:tspan % at each time point
    
    for nn = 1:2^nbd % each profile contains 2^n possible configurations, for each configurati
        
        cur_num_arr = zeros(nbd,1);
        cur_conc_arr = zeros(nbd,1);
            
        for bb = 1:nbd % each configuration has 3 binding boxes
            
            if bins(nn,bb) == 0  % if nothing binds
                cur_num_arr(bb) = 0; % # of TF = 0
                cur_conc_arr(bb) = 0; % conc of protein doesn't matter
            end
            
            if bins(nn,bb) == 1 % if this binding box is bound
                cur_num_arr(bb) = 1; % # of TF = 1 (for now)
                cur_conc_arr(bb) = tfunc(tt); % conc of protein linearly vary with time
            end
        
        end
        
        % record two fields for each configuration to complete information
        % for one configuration
        cur_config(nn).num_arr = cur_num_arr;
        cur_config(nn).conc_arr = cur_conc_arr;
        
    end  
    
    t_overall = [t_overall,cur_config];
    
end



                
    