load_data;

time=2:10;


t= p_value("ps_"+"123");err=p_std("ps_"+"123");
s=errorbar(time,t,err,'ko-',"LineWidth",2); hold on;

t= p_value("ps_"+"1");err=p_std("ps_"+"1");
s=errorbar(time,t,err,'ro-',"LineWidth",2); hold on;
t= p_value("ps_"+"2");err= p_std("ps_"+"2");
s=errorbar(time,t,err,'go-',"LineWidth",2); hold on;
t= p_value("ps_"+"3");err=p_std("ps_"+"3");
s=errorbar(time,t,err,'bo-',"LineWidth",2); hold on;
t= p_value("ps_");err=p_std("ps_");
s=errorbar(time,t,err,'mo-',"LineWidth",2); hold on;
legend('ps','1*2*','1*3*','2*3*','1*2*3*')
%%
figure;


t= p_value("ps_"+"123");err=p_std("ps_"+"123");
s=errorbar(time,t,err,'ko-',"LineWidth",2); hold on;

t= p_value("ps_"+"12");err=p_std("ps_"+"12");
s=errorbar(time,t,err,'ro-',"LineWidth",2); hold on;
t= p_value("ps_"+"13");err= p_std("ps_"+"13");
s=errorbar(time,t,err,'go-',"LineWidth",2); hold on;
t= p_value("ps_"+"23");err=p_std("ps_"+"23");
s=errorbar(time,t,err,'bo-',"LineWidth",2); hold on;

legend('ps','1*','2*','3*')


%%
figure;

t= p_value("pv_"+"123");err=p_std("pv_"+"123");
s=errorbar(time,t,err,'ko-',"LineWidth",2); hold on;

t= p_value("pv_"+"1");err=p_std("pv_"+"1");
s=errorbar(time,t,err,'ro-',"LineWidth",2); hold on;
t= p_value("pv_"+"2");err= p_std("pv_"+"2");
s=errorbar(time,t,err,'go-',"LineWidth",2); hold on;
t= p_value("pv_"+"3");err=p_std("pv_"+"3");
s=errorbar(time,t,err,'bo-',"LineWidth",2); hold on;
t= p_value("pv_");err=p_std("pv_");
s=errorbar(time,t,err,'mo-',"LineWidth",2); hold on;
legend('pv','1*2*','1*3*','2*3*','1*2*3*')
%%
figure;
t= p_value("pv_"+"123");err=p_std("pv_"+"123");
s=errorbar(time,t,err,'ko-',"LineWidth",2); hold on;

t= p_value("pv_"+"12");err=p_std("pv_"+"12");
s=errorbar(time,t,err,'ro-',"LineWidth",2); hold on;
t= p_value("pv_"+"13");err= p_std("pv_"+"13");
s=errorbar(time,t,err,'go-',"LineWidth",2); hold on;
t= p_value("pv_"+"23");err=p_std("pv_"+"23");
s=errorbar(time,t,err,'bo-',"LineWidth",2); hold on;

legend('pv','1*','2*','3*')