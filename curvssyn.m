
function [t,v,w,sq]=curvssyn(dt,gbar,T,Tfin)

VCl = -68;	    % mV
A = 4*pi*1e-6;	% cm^2 patch area
Cm = 1;         % micro F/cm^2
gCl = 0.3;      % mS/cm^2
tau = Cm/gCl;   % ms
E = 60;
z = 2*Cm/dt;

%Tfin = 10*T;
Nt = ceil(1+Tfin/dt);                  % number of time steps
v = zeros(Nt,1);  w = v; t = v;  sq = v; % preallocate space

j = 1;
g0 = 0;

for j=2:Nt,

   t(j) = (j-1)*dt;

   sq(j) = (mod(t(j),T)<1); 
   g1 = gbar*sq(j);
   
   v(j) = ( (z-gCl-g0)*v(j-1) + (g0+g1)*E ) / (z+gCl+g1);
   
   w(j) = ( (z-gCl)*w(j-1) + (g0+g1)*E )/ (z+gCl);
   
   g0 = g1;

end

