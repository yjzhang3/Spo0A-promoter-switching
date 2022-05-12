function [n_arr,c_arr] = gen_config_prof(config,TF_conc,RNAp_conc)
% input: a specific configuration (specified by 01 array)
% output: the number and concentration profile for each configuration

n_arr = zeros(2,1); % assume only two types of protein, TF and RNAP
c_arr = zeros(2,1); % assume only two types of protein, TF and RNAP

n_arr(end) = sum(config(end-1:end)); % number of rnap is the total number
% of bound promoters (site 4th and 5th)
c_arr(end) = RNAp_conc; 

% assume besides RNAP, all the other binding sites are bound by *one* type
% of TF
n_arr(1) = sum(config(1:3));
c_arr(1) = TF_conc;

end
    