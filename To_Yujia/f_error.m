function err = f_error(p_data,args)
t=num2cell(args);

[   k1,k2,k3,k12,k13,k23,k123, ...
    k0a,k1a,k2a,k3a,k12a,k13a,k23a,k123a, ...
    k0h,k1h,k2h,k3h,k12h,k13h,k23h,k123h,...
    A,H,kscale] = deal(t{:});

t_args=[
    k1,k2,k3,k12,k13,k23,k123, ...
    k0a,k1a,k2a,k3a,k12a,k13a,k23a,k123a, ...
    k0h,k1h,k2h,k3h,k12h,k13h,k23h,k123h];

strains=["","1","2","3","12","13","23","123"];
muts={[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]};
muts_map= containers.Map(strains,muts);

aps=get_aps(2:10);

err=0;
for strain = strains
    v="pv_"+strain;
    pv=f_dynamics(muts_map(strain),t_args,'v',aps,A);
    err=err+immse(pv'.*kscale,p_data(v));
end

for strain = strains
    s="ps_"+strain;
    ps=f_dynamics(muts_map(strain),t_args,'s',aps,H);
    err=err+immse(ps'.*kscale,p_data(s));
end

end

