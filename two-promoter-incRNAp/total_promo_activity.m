load('A_conc_t_reasonable.mat', 'A_conc_t')
load('H_conc_t_reasonable.mat', 'H_conc_t')

nbd = 5;
energyi = [5,8,-9,0,0,-1.8,1,3,5,2.4,4,-8,6.2,7,3.5];
TF_conc_t = get_aps(2:10);
RNAp_conc_t_Ps = H_conc_t;
RNAp_conc_t_Pv = A_conc_t;
mut_mat = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[0,0,1];[0,0,0]];
vmax = 450;

%%
TR_overall = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNAp_conc_t_Ps,RNAp_conc_t_Pv,mut_mat(3,:),vmax);