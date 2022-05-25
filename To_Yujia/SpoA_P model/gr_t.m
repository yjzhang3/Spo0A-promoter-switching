function gr_t = gr_t(t)

para=[       10.1572    1.95 120.5862    1.6943    0.622   13.2283    3.0000    0.3000 ];
kc=para(1);
vm=para(2);
Kn=para(3);
h=para(4);kd=para(5);Kd=para(6);m=para(7);fn=para(8);

C0=0.1;

[tts,y]=ode15s(@(t,nc)F(t,nc,para(1:8)),[0,30],[100,C0]');
%  plot(tts,y(:,1),'k','Linewidth',2); hold on;
n=y(:,1);

dc=vm.*(n.^h./(n.^h+Kn.^h));
gr_t=interp1(tts,dc,t);

end
