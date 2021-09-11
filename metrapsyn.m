function [t,V,g] = metrapsyn(dt,Tfin,gsyn)

VCl = -68;	    % mV
Cm = 1;         % micro F/cm^2
gCl = 0.3;      % mS/cm^2
z = 2/dt;

Nt = ceil(1+Tfin/dt);                    % number of time steps
V = zeros(Nt,1);  t = V;  	 % preallocate space
g = zeros(Nt,length(gsyn.gmax));

j = 1;
V(1) = VCl;
a0 = gCl/Cm;
b0 = gCl*VCl/Cm;

for j=2:Nt

   t(j) = (j-1)*dt;

   g(j,:) = gsyn.gmax.*((t(j)-gsyn.t1)./gsyn.taua).*exp(1-(t(j)-gsyn.t1)./gsyn.taua).*(t(j)>gsyn.t1);
   a1 = (gCl+sum(g(j,:)))/Cm;
   tmp = g(j,:)*(gsyn.Vsyn)';
   b1 = (gCl*VCl + g(j,:)*(gsyn.Vsyn)')/Cm;
   
   V(j) = ( (z-a0)*V(j-1) + b0+b1 ) / (z+a1);
   
   a0 = a1;
   b0 = b1;

end

