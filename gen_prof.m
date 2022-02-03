function [n_arr,c_arr] = gen_prof(config,TF_conc)
% input: a specific configuration (specified by 01 array)
% output: the number and concentration profile for each configuration

n_arr = zeros(length(config)*10,1); % initialize an array as large as possible
c_arr = zeros(length(config)*10,1); % initialize an array as large as possible

if config(end) == 1 % assume RNAP conc is always 0.1 and of one molecule
    n_arr(end) = 1;
    c_arr(end) = 1; 
end

% assume besides RNAP, all the other binding sites are bound by *one* type
% of TF
n_arr(1) = sum(config(1:end-1));
c_arr(1) = TF_conc;

end
    