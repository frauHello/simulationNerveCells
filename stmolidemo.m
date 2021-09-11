

function [t,V,INa,IK,IA]=stmolidemo(dt,tfin,I0)

Nt = tfin/dt + 1;
V = zeros(1,Nt);t = V;Istim = V;INa = V;IK = V;IA = V;
A = 4*pi*(1e-6);       % cm^2
Cm = 1.5;
g = struct('K',7,'Na',30,'L',1,'A',16);   
E = struct('K',-90,'Na',45,'L', -70);  
Vr = fsolve(@(v) irest(v,g,E),-120);

V(1) = Vr;
h = hinf(Vr);
n = ninf(Vr);
b = binf(Vr);

for j=2:Nt

    tmp = 2*tauh(V(j-1)); 
    h = ( h*(tmp-dt) + 2*dt*hinf(V(j-1)) ) / ( dt + tmp );

    tmp = 2*taun(V(j-1)); 
    n = ( n*(tmp-dt) + 2*dt*ninf(V(j-1)) ) / ( dt + tmp );

    tmp = 2*taub(V(j-1)); 
    b = ( b*(tmp-dt) + 2*dt*binf(V(j-1)) ) / ( dt + tmp );

    a = ainf(V(j-1));
    m = minf(V(j-1));
    
    t(j) = dt*(j-1);
    
    Istim = (t(j)-dt/2>5)*(t(j)-dt/2<7)*I0*1e-6;
    
    top = 2*Cm*V(j-1)/dt + g.Na*m*h*E.Na + g.K*n*E.K + g.A*a*b*E.K + g.L*E.L + Istim/A;
    bot = 2*Cm/dt + g.Na*m*h + g.K*n + g.A*a*b + g.L;

    Vmid = top/bot;

    V(j) = 2*Vmid - V(j-1);
    
    INa(j) = g.Na*minf(V(j))*h*(V(j)-E.Na);
    IK(j) = g.K*n*(V(j)-E.K);
    IA(j) = g.A*ainf(V(j))*b*(V(j)-E.K);

end

%subplot(1,2,1)
%plot(t,V,'k')
%xlim([0 15])
%xlabel('t  (ms)','fontsize',14)
%ylabel('V  (mV)','fontsize',14)
%subplot(1,2,2)
%plot(t,IA,'k')
%hold on
%plot(t,IK,'r--')
%ylabel('\muA/cm^2','fontsize',14)
%plot(t,INa,'r')
%xlim([6 9])
%ylabel('\muA/cm^2','fontsize',14)
%legend('I_A','I_K','I_{Na}')
%hold off
%xlabel('t  (ms)','fontsize',14)
return

function val = irest(v,g,E)
val = g.Na*minf(v)*hinf(v)*(v-E.Na) + g.K*ninf(v)*(v-E.K) + ...
      g.A*ainf(v)*binf(v)*(v-E.K) + g.L*(v-E.L);

% Na functionals

function val = minf(v)
val = 1./(1+exp(-(v+35)/4));

function val = hinf(v)
val = 1./(1+exp((v+35)/4));

function val = tauh(v)
val = 2*232*28./(4*pi*(v+74).^2+28^2) - 0.15;

% K functionals

function val = ninf(v)
val = 1./(1+exp(-(v+35)/4));

function val = taun(v)
val = 0.5;

% A type K functionals

function val = ainf(v)
val = 1./(1+exp(-(v+27)/8.8));

function val = binf(v)
val = 1./(1+exp((v+68)/6.6));

function val = taub(v)
val = 15;

