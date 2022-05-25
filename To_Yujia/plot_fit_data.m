load_data;
load("full_sols.mat");
time=2:10;


args=full_sols{2}{1};
t=num2cell(args);

[   k1,k2,k3,k12,k13,k23,k123, ...
    k0a,k1a,k2a,k3a,k12a,k13a,k23a,k123a, ...
    k0h,k1h,k2h,k3h,k12h,k13h,k23h,k123h,...
    A,H,kscale] = deal(t{:});

t_args=[
    k1,k2,k3,k12,k13,k23,k123, ...
    k0a,k1a,k2a,k3a,k12a,k13a,k23a,k123a, ...
    k0h,k1h,k2h,k3h,k12h,k13h,k23h,k123h];


t= p_value("ps_"+"123");err=p_std("ps_"+"123");
s=errorbar(time,t,err,'ko'); hold on;

t= p_value("ps_"+"1");err=p_std("ps_"+"1");
s=errorbar(time,t,err,'ro'); hold on;
t= p_value("ps_"+"2");err= p_std("ps_"+"2");
s=errorbar(time,t,err,'go'); hold on;
t= p_value("ps_"+"3");err=p_std("ps_"+"3");
s=errorbar(time,t,err,'bo'); hold on;
t= p_value("ps_");err=p_std("ps_");
s=errorbar(time,t,err,'mo'); hold on;
legend('ps','1*2*','1*3*','2*3*','1*2*3*')

aps=get_aps(2:10);
s1=f_dynamics([1,1,1],t_args,'s',aps,H);
s2=f_dynamics([1,0,0],t_args,'s',aps,H);
s3=f_dynamics([0,1,0],t_args,'s',aps,H);
s4=f_dynamics([0,0,1],t_args,'s',aps,H);
s5=f_dynamics([0,0,0],t_args,'s',aps,H);
plot(time,kscale*s1,'k'); hold on;
plot(time,kscale*s2,'r'); hold on;
plot(time,kscale*s3,'g'); hold on;
plot(time,kscale*s4,'b'); hold on;
plot(time,kscale*s5,'m'); hold on;
%%
figure;


t= p_value("ps_"+"123");err=p_std("ps_"+"123");
s=errorbar(time,t,err,'ko'); hold on;

t= p_value("ps_"+"12");err=p_std("ps_"+"12");
s=errorbar(time,t,err,'ro'); hold on;
t= p_value("ps_"+"13");err= p_std("ps_"+"13");
s=errorbar(time,t,err,'go'); hold on;
t= p_value("ps_"+"23");err=p_std("ps_"+"23");
s=errorbar(time,t,err,'bo'); hold on;

legend('ps','1*','2*','3*')
s1=f_dynamics([1,1,1],t_args,'s',aps,H);
s2=f_dynamics([1,1,0],t_args,'s',aps,H);
s3=f_dynamics([1,0,1],t_args,'s',aps,H);
s4=f_dynamics([0,1,1],t_args,'s',aps,H);
plot(time,kscale*s1,'k'); hold on;
plot(time,kscale*s2,'r'); hold on;
plot(time,kscale*s3,'g'); hold on;
plot(time,kscale*s4,'b'); hold on;

%%
figure;

t= p_value("pv_"+"123");err=p_std("pv_"+"123");
s=errorbar(time,t,err,'ko'); hold on;
t= p_value("pv_"+"1");err=p_std("pv_"+"1");
s=errorbar(time,t,err,'ro'); hold on;
t= p_value("pv_"+"2");err= p_std("pv_"+"2");
s=errorbar(time,t,err,'go'); hold on;
t= p_value("pv_"+"3");err=p_std("pv_"+"3");
s=errorbar(time,t,err,'bo'); hold on;
t= p_value("pv_");err=p_std("pv_");
s=errorbar(time,t,err,'mo'); hold on;
legend('pv','1*2*','1*3*','2*3*','1*2*3*')

s1=f_dynamics([1,1,1],t_args,'v',aps,A);
s2=f_dynamics([1,0,0],t_args,'v',aps,A);
s3=f_dynamics([0,1,0],t_args,'v',aps,A);
s4=f_dynamics([0,0,1],t_args,'v',aps,A);
s5=f_dynamics([0,0,0],t_args,'v',aps,A);

plot(time,kscale*s1,'k'); hold on;
plot(time,kscale*s2,'r'); hold on;
plot(time,kscale*s3,'g'); hold on;
plot(time,kscale*s4,'b'); hold on;
plot(time,kscale*s5,'m'); hold on;


%%
figure;
t= p_value("pv_"+"123");err=p_std("pv_"+"123");
s=errorbar(time,t,err,'ko'); hold on;

t= p_value("pv_"+"12");err=p_std("pv_"+"12");
s=errorbar(time,t,err,'ro'); hold on;
t= p_value("pv_"+"13");err= p_std("pv_"+"13");
s=errorbar(time,t,err,'go'); hold on;
t= p_value("pv_"+"23");err=p_std("pv_"+"23");
s=errorbar(time,t,err,'bo'); hold on;

legend('pv','1*','2*','3*')

s1=f_dynamics([1,1,1],t_args,'v',aps,A);
s2=f_dynamics([1,1,0],t_args,'v',aps,A);
s3=f_dynamics([1,0,1],t_args,'v',aps,A);
s4=f_dynamics([0,1,1],t_args,'v',aps,A);

plot(time,kscale*s1,'k'); hold on;
plot(time,kscale*s2,'r'); hold on;
plot(time,kscale*s3,'g'); hold on;
plot(time,kscale*s4,'b'); hold on;


