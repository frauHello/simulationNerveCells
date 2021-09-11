
stim = struct('amp',5e-6,'t1',1,'t2',21,'Tfin',40);
[t,V] = stE(0.01,stim);
plot(t,V,'k')
box off
xlabel('t  (ms)','fontsize',14)
ylabel('V  (mV)','fontsize',14)

