function prom = total_promo_activity(config,a,b,c)

nbd = 4;
energyi = ?
TF_conc_t = get_aps(2:10);
RNAp_conc_t_Ps = ?
RNAp_conc_t_Pv = ?
mut_mat = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[0,0,1];[0,0,0]];
vmax = ?

Ps_TR = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNAp_conc_t,mut,vmax);
Pv_TR = time_dep_TR_new_wSigma(nbd,energyi,TF_conc_t,RNAp_conc_t,mut,vmax);

prom = a*Ps_TR*config(1)+b*Pv_TR*config(2)+c*Ps_TR*Pv_TR*config(1)*config(2);