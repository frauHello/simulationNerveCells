

function [t,V] = stE2d(dt,Tfin,I0)

A = 4*pi*(1e-6);% (cm)^2

E = struct('K', -77, 'Na', 56, 'Cl', -68); % reversal potentials, mV

G = struct('K', 36, 'Na', 120, 'Cl', 0.3);  % channel conductances, mS/cm^2

Vr = fsolve(@(V) Iss(V,E,G),-71)   % find rest potential 

Nt = ceil(Tfin/dt);
t = zeros(Nt,1); V = t; n = V;
t(1) = 0;
V(1) = Vr;
n(1) = an(Vr)/(an(Vr)+bn(Vr));  

td = 1/dt;

for j = 2:Nt;

      t(j) = (j-1)*dt;

      Istim = I0*(t(j)-dt/2>2)*(1e-6);
      
      a = an(V(j-1));  b = bn(V(j-1)); c = (a+b)/2;
      n(j) = ( (td-c)*n(j-1) + a ) / (td + c);
      
      cK = G.K*n(j)^4;
      
      a = am(V(j-1));  b = bm(V(j-1)); 
      m = a / (a + b);
      
      h = 0.87 - n(j);
      %h = 0.7 - n(j)^2;
      
      cNa = G.Na*m^3*h;
      
      top = 2*V(j-1)*td + cK*E.K + cNa*E.Na + G.Cl*E.Cl + Istim/A;
      
      bot = 2*td + cK + cNa + G.Cl;
      
      Vmid = top/bot;

      V(j) = 2*Vmid - V(j-1); 

end

%subplot(2,1,1)
%plot(t,V)
%ylim([-82 60])
%box off
%xlabel('t  (ms)','fontsize',14)
%ylabel('V  (mV)','fontsize',14)

%subplot(2,1,2)
%plot(V,n)
%ylim([.25 .85])
%box off
%xlabel('V  (mV)','fontsize',14)
%ylabel('n','fontsize',14)

function val = Iss(V,E,G)
m = am(V)/(am(V)+bm(V));
n = an(V)/(an(V)+bn(V));
%h = 0.87 - n;
h = 0.7 - n.^2;
val = G.Na*m^3*h*(V-E.Na) + G.K*n^4*(V-E.K) + G.Cl*(V-E.Cl);

function val = an(V)
val = .01*(10-(V+71))./(exp(1-(V+71)/10)-1);

function val = bn(V)
val = .125*exp(-(V+71)/80);

function val = am(V)
val = .1*(25-(V+71))./(exp(2.5-(V+71)/10)-1);

function val = bm(V)
val = 4*exp(-(V+71)/18);

function val = ah(V)
val = 0.07*exp(-(V+71)/20);

function val = bh(V)
val = 1./(exp(3-(V+71)/10)+1);
