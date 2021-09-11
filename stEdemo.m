
function [t,V,n,m,h,INa,IK,ICl] = stEdemo(dt,Tfin,I0)

A = 4*pi*(1e-6);% (cm)^2

E = struct('K', -77, 'Na', 56, 'Cl', -68); % reversal potentials, mV

G = struct('K', 36, 'Na', 120, 'Cl', 0.3);  % channel conductances, mS/cm^2

Vr = fsolve(@(V) Iss(V,E,G),-71);   % find rest potential 

Nt = ceil(Tfin/dt);
t = zeros(Nt,1); V = t; n = t; m = t; h = t; INa = t; IK = t; ICl = t;
t(1) = 0;
V(1) = Vr;
n(1) = an(Vr)/(an(Vr)+bn(Vr));  
m(1) = am(Vr)/(am(Vr)+bm(Vr));  
h(1) = ah(Vr)/(ah(Vr)+bh(Vr));  

td = 1/dt;

for j = 2:Nt;

      t(j) = (j-1)*dt;

      Istim = I0*(t(j)-dt/2>2)*(t(j)-dt/2<4)*(1e-6);
      
      a = an(V(j-1));  b = bn(V(j-1)); c = (a+b)/2;
      n(j) = ( (td-c)*n(j-1) + a ) / (td + c);
      
      cK = G.K*n(j)^4;
      
      a = am(V(j-1));  b = bm(V(j-1)); c = (a+b)/2;
      m(j) = ( (td-c)*m(j-1) + a ) / (td + c);
      
      a = ah(V(j-1));  b = bh(V(j-1)); c = (a+b)/2;
      h(j) = ( (td-c)*h(j-1) + a ) / (td + c);
      
      cNa = G.Na*m(j)^3*h(j);
      
      top = 2*V(j-1)*td + cK*E.K + cNa*E.Na + G.Cl*E.Cl + Istim/A;
      
      bot = 2*td + cK + cNa + G.Cl;
      
      Vmid = top/bot;

      V(j) = 2*Vmid - V(j-1); 
      
      ICl(j) = G.Cl*(V(j)-E.Cl);

      IK(j) = G.K*n(j)^4*(V(j)-E.K);

      INa(j) = G.Na*m(j)^3*h(j)*(V(j)-E.Na);

end


function val = Iss(V,E,G)
val = G.Na*(am(V)/(am(V)+bm(V)))^3*(ah(V)/(ah(V)+bh(V)))*(V-E.Na) + ...
      G.K*(an(V)/(an(V)+bn(V)))^4*(V-E.K) + G.Cl*(V-E.Cl);

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
