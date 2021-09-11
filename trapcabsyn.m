
function trapcabsyn(cab,stim,pinc)

Cm = 1;		% micro F / cm^2
g_L = 1/15; %0.3;     		% mS / cm^2
R_2 = 0.3; % 0.034;		% k Ohm cm
dx = cab.dx;
dt = cab.dt;
Nx = cab.ell/dx;                % patch length
A = 2*pi*cab.rad*dx;            % patch surface area
x = dx/2:dx:cab.ell-dx/2;       % vector of patch midpoints

Nt = ceil(stim.Tfin/dt);

v = zeros(Nx,1);		% initial conditions

e = ones(Nx,1);
S = spdiags([-e 2*e -e], -1:1, Nx, Nx)/dx/dx;
S(1,1) = 1/dx/dx;
S(Nx,Nx) = 1/dx/dx;

tau = Cm/g_L;
lambda = sqrt(cab.rad/(2*R_2*g_L));
A = 2*pi*cab.rad*dx;

e1 = zeros(Nx,1);
eloc = round(Nx*stim.loc/cab.ell);
vsyn = 70;
Iapp = stim.Gsyn*(dt/2)/A/Cm;

B = (dt/2)*(speye(Nx)+lambda^2*S)/tau;
B = speye(Nx) + B;
bloc = (eloc-1)*(Nx+1)+1;
dBe = B(bloc);
figure('Name',' Response of the passive cable to ?-function synaptic input with gsyn=100 nS','Position',[1 1 1400 770]);

t = 0;
curve=animatedline();
addpoints(curve,x,t*e,v);

head=scatter3(x,t*e,v,'filled');
drawnow;
pause(0.01);
delete(head);

hold on

c0 = Iapp.*(t./stim.tau).*exp(1-t./stim.tau).*(t>stim.t1);
t = dt;
c1 = Iapp.*(t./stim.tau).*exp(1-t./stim.tau).*(t>stim.t1);

r = zeros(Nx,1);
r(eloc) = vsyn*(c0 + c1)';

vhot = zeros(length(eloc),Nt);
cc=jet(Nt);
for j=2:Nt,
    
    B(bloc) = dBe + c1;

    v = B\r; 

    vhot(:,j) = v(eloc);

    if mod(j,pinc) == 0
        curve=animatedline('color',cc(j,:));
  addpoints(curve,x,t*e,v);
  head=scatter3(x,t*e,v,'filled');
  drawnow;
  pause(0.01);
  delete(head);
      
    end

    t = j*dt;
    c0 = c1;
    c1 = Iapp.*((t-stim.t1)./stim.tau).*exp(1-(t-stim.t1)./stim.tau).*(t>stim.t1);

    r = 2*v - r;
    r(eloc) = r(eloc) + vsyn*(c0 + c1)';

end

xlabel('x (cm)','fontsize',14)
ylabel('t (ms)','fontsize',14)
zlabel('v (mV)','fontsize',14)

hold off


