
function [IC,ICl,v,t]=bepswI(dt,Tfin,I,time)

VCl = -68;	    % mV
A = 4*pi*1e-6;	% cm^2 patch area
Cm = 1;         % micro F/cm^2
gCl = 0.3;      % mS/cm^2
tau = Cm/gCl;   % ms

Nt = ceil(1+Tfin/dt);                  % number of time steps
v = zeros(Nt,1);  t = v; Istim = v;    % preallocate space
IC = v;
ICl = v;

v(1) = VCl;                 % initialize v

for j=2:Nt

   t(j) = j*dt;

   Istim(j) = (t(j)>2)*(t(j)<(2+time))*1e-6*I;    % 10 pA 20 ms pulse

   v(j) = (v(j-1) + dt*(VCl/tau + Istim(j)/A/Cm))/(1+dt/tau);  % back euler

   ICl(j) = gCl*(v(j)-VCl);

   IC(j) = Istim(j)/A - ICl(j);

end

