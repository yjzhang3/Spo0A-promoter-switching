function s = fss_fit(no,gr,va,vc,vs,kg,kt,nt,p0e)
%no: species order
%out:steady-state value
n=config('Spo0A_fit.txt');

xi=zeros(1,n);
 [~,y]=ode15s(@(t,x)Spo0A_full(t,x,gr,va,vc,vs,kg,kt,nt,p0e),[0,80],xi);
 s=[];
 for i=1:length(no)
s=[s y(end,no(i))];
 end
end
