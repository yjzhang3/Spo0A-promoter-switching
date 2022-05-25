function dnc=F(t,nc,para)
kc=para(1);
vm=para(2);
Kn=para(3);
h=para(4);kd=para(5);Kd=para(6);m=para(7);fn=para(8);
% n-nutrient   c-cell    dnc-[dn,dc] nc-[n,c]
%nct-nc before the delay
n=nc(1);
c=nc(2);

dn=-c.*kc.*(vm.*(n.^h./(n.^h+Kn.^h)))+c.*(kd*Kd.^m./(Kd^m+n.^m))*kc*fn;
dc=c.*(vm*(n.^h./(n.^h+Kn.^h))-kd*Kd.^m./(Kd^m+n.^m));
% kc:the nutrient consuming rate
% vm: maximum growth rate
% Kn: half-maximum level of nutrient for growth
% kd: death rate
% Kd/m : suppose the death rate is determained by nutrient level in hill
% manner
% kn2 k2 K2: suppose there is a negative effect of the cell density on the
% generic growth rate
dnc=[dn dc]';
end