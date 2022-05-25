function diff = objective_function_together(nbd,p,TF_conc_t,H_conc_t,A_conc_t,mut_mats,mut_matv,real_datas,real_datav,n_strains,n_strainv)
% this objective function is for validation or real fitting when time-dependent RNAp is
% known

% n_strain specifies the number of strains used for fitting

% assign parameters (RNAp is not independent of energies, so it's not a
% parameter to optimize)
% 1,2,3,s,v,12,13,23,1s,2s,3s,1v,2v,3v
single_energyi = p(1:5);
TF_int_energyi = p(6:8);
Ps_int_energyi = p(9:11);
Pv_int_energyi = p(12:14);
vmax = p(15:end);

% 1,2,3,p,12,13,1p,23,2p,3p
energyi_s = zeros(10,1);
energyi_s(1:4) = single_energyi(1:4);
energyi_s([5,6,8]) = TF_int_energyi;
energyi_s([7,9,10]) = Ps_int_energyi;

energyi_v = zeros(10,1);
energyi_v(1:4) = single_energyi([1:3,5]);
energyi_v([5,6,8]) = TF_int_energyi;
energyi_v([7,9,10]) = Pv_int_energyi;

vmax_s = vmax(1:(2^(nbd-1)));
vmax_v = vmax((2^(nbd-1)+1):(2*2^(nbd-1)));

% all data follow this format:
% each row: every hour
% each column: every mutant

% generate simulated data for all WT and mutated type

    
parfor mm = 1:n_strains
    sim_data_s(:,mm) = time_dep_TR_new_wSigma(nbd,energyi_s,TF_conc_t,H_conc_t,mut_mats(mm,:),vmax_s);
end

parfor nn = 1:n_strainv
    sim_data_v(:,nn) = time_dep_TR_new_wSigma(nbd,energyi_v,TF_conc_t,A_conc_t,mut_matv(nn,:),vmax_v);
end

sim_data = [sim_data_s;sim_data_v];
real_data = [real_datas;real_datav];

% now sim_data should be the same dimension as real_data

diff = weighted_msd(real_data(:),sim_data(:));

end