
function [t,V]=ml2(dt, Tfin, IsA)

Cm = 1;
g = struct('K',20,'Ca',15,'Cl',5,'syn',30);
E = struct('K',-80,'Ca',100,'Cl',-50,'syn',-80);

Nt = ceil(Tfin/dt);
t = linspace(0,Tfin,Nt);
V = zeros(Nt,2); w = V;

V(1,:) = [-10 -20];
w(1,:) = minf(V(1,:));

for j=2:Nt,

    tmpt = tauw(V(j-1,:));
    tmpm = minf(V(j-1,:));
    tmps = sinf(fliplr(V(j-1,:)));
  
    w(j,:) = (tmpt.*w(j-1,:) + tmpm*dt)./(dt + tmpt);

    top = Cm*V(j-1,:) + dt*(g.Ca*tmpm*E.Ca + g.K*w(j,:)*E.K + ...
                            g.Cl*E.Cl + g.syn*tmps*E.syn + IsA*[1 1]);
    bot = Cm + dt*(g.Ca*tmpm + g.K*w(j,:) + g.Cl + g.syn*tmps);

    V(j,:) = top./bot;

end












function val = minf(V)
val = (1 + tanh(V/15))/2;
return

function val = tauw(V)
val = 1./cosh(V/30);
return

function val = sinf(V)
val = (1 + tanh(V/15))/2;
return
