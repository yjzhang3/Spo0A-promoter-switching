nbd = 3;
bins = dec2bin(0:(2^nbd-1), nbd) - '0';
tspan = 10;

t_overall = [];
for tt = 1:tspan % at each time point

    for nn = 1:2^nbd % each profile contains 2^n possible configurations, for each configurati

        cur_num_arr = zeros(nbd,1);
        cur_conc_arr = zeros(nbd,1);

        for bb = 1:nbd % each configuration has 3 binding boxes

            if bins(nn,bb) == 0  % if nothing binds
                cur_num_arr(bb) = 0; % # of TF = 0
                cur_conc_arr(bb) = 0; % conc of protein doesn't matter
            else % if this binding box is bound
                cur_num_arr(bb) = 1; % # of TF = 1 (for now)
                cur_conc_arr(bb) = 10; % conc of protein linearly vary with time
            end

        end

        cur_config(nn).num_arr = cur_num_arr;
        cur_config(nn).conc_arr = cur_conc_arr;
    end 
    t_overall = [t_overall,cur_config];
end