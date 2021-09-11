

function [t,gE1,V1,gE2,V2]=twocell(dt,Tfin,P)

Nt = ceil(Tfin/dt);

tauE = 2; VE = 0;

gL = 0.3; VL = -68; Cm = 1;

winp = 0.5; w21 = 0.5;

aE = (2*tauE-dt)/(2*tauE+dt);  bE = 2/(2*tauE+dt);

tref = 3;
Vthr = -50;
Vres = -70;

gE1 = zeros(Nt,1); gE2 = gE1; 
V1 = VL*ones(Nt,1); V2 = V1;

sp1 = 0; sp2 = 0; T1 = -10; T2 = -10;

t = gE1; 

for j=1:Nt-1,

    t(j+1) = j*dt;

    spin = (t(j+1)/P==round(t(j+1)/P));

    gE1(j+1) = aE*gE1(j) + bE*winp*spin;

    gE2(j+1) = aE*gE2(j) + bE*w21*sp1;

    if t(j+1)-T1 > tref,

       n1 = (2*Cm/dt-(gL+gE1(j)))*V1(j) + 2*gL*VL;
       d1 = 2*Cm/dt + gL  + gE1(j+1);
       tmp = n1/d1;

       if tmp < Vthr
          V1(j+1) = tmp;
       else
          V1(j+1) = Vres;
          sp1 = 1;
          T1 = t(j+1);
       end

    else

       sp1 = 0;
       V1(j+1) = Vres;

    end

    if t(j+1)-T2 > tref,

       n2 = (2*Cm/dt-(gL+gE2(j)))*V2(j) + 2*gL*VL;
       d2 = 2*Cm/dt + gL +gE2(j+1);
       tmp = n2/d2;

       if tmp < Vthr
          V2(j+1) = tmp;
       else
          V2(j+1) = Vres;
          sp2 = 1;
          T2 = t(j+1);
       end

    else

       sp2 = 0;
       V2(j+1) = Vres;

    end
    
end  % for j



