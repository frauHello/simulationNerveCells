
function stEthresh(Ilo,Ihi,Iinc)

stim = struct('t1',1,'t2',3,'amp',1e-6,'Tfin',10);

I0 = Ilo:Iinc:Ihi;

for j = 1:length(I0)

    stim.amp = I0(j);
    [t,V] = stE(0.01,stim);
    [vmax(j),tj] = max(V);
    tmax(j) = t(tj);

end


[ax,h1,h2] = plotyy(I0*10^6,vmax,I0*10^6,tmax);
set(ax(1),'ycolor','k','ylim',[-80 60],...
    'ytick',[-80 -60 -40 -20 0 20 40 60])
set(ax(2),'ycolor','r','ylim',[0 10],...
    'ytick',[0 2 4 6 8 10])
set(h1,'color','k')
set(h2,'color','r')
xlabel('I_{0} (pA)','fontsize',14)
ylabel(ax(1),'V_{max} (mV)','fontsize',14)
ylabel(ax(2),'t_{max} (ms)','fontsize',14)
box off

return
