

function [t,gE,V]=threecell(dt,Tfin,P)

Nt = ceil(Tfin/dt);

tauE = 2; VE = 0;

gL = 0.3; VL = -68; Cm = 1;

Win = [0.5 0 0]';

w = 0.5;   %0.26;
W = [0 0 0; 0.5 0 0; w w 0];

aE = (2*tauE-dt)/(2*tauE+dt);  bE = 2/(2*tauE+dt);

tref = 3;
Vthr = -50;
Vres = -70;

gE = zeros(3,Nt); 
V = VL*ones(3,Nt);
t = zeros(Nt,1);

sp = zeros(3,1); T = -10*ones(3,1); ref = zeros(3,1);

for j=1:Nt-1,

    t(j+1) = j*dt;

    spin = (t(j+1)/P==round(t(j+1)/P));

    gE(:,j+1) = aE*gE(:,j) + bE*(W*sp + spin*Win);

    top = (2*Cm/dt - (gL+gE(:,j))).*V(:,j) + 2*gL*VL;
    bot = 2*Cm/dt + gL + gE(:,j+1);
    V(:,j+1) = top./bot;

    ref = find( t(j+1)-T < tref );
    V(ref,j+1) = Vres;

    sp = (V(:,j+1)>Vthr);

    if sum(sp)>0
        fsp = find(sp>0);
        T(fsp) = t(j+1);
        V(fsp,j+1) = Vres;
    end


end  % for j

