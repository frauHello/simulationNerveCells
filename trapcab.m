

function trapcab(cab,stim,pinc)

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

B = (speye(Nx)+lambda^2*S)/tau;

[L,U] = lu(speye(Nx)+B*dt/2);

e1 = zeros(Nx,1);
eloc = round(Nx*stim.loc/cab.ell);
Iapp = stim.amp*(dt/2)/A/Cm;

vhot = zeros(length(eloc),Nt);
figure('Name','Response of the single passive cable to  double current injection','Position',[1 1 1400 770]);
t = 0;
curve=animatedline();
addpoints(curve,x,t*e,v);

head=scatter3(x,t*e,v,'filled');
drawnow;
pause(0.01);
delete(head);
hold on

f0 = Iapp.*(t>stim.t1).*(t<stim.t2);
t = dt;
f1 = Iapp.*(t>stim.t1).*(t<stim.t2);

r = zeros(Nx,1);
r(eloc) = (f0 + f1)';
cc=jet(Nt);
for j=2:Nt,
    
    v = U \ ( L \ r );

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
    f0 = f1;
    f1 = Iapp.*(t>stim.t1).*(t<stim.t2);

    r = 2*v - r;
    r(eloc) = r(eloc) + (f0 + f1)';

end
xlabel('x (cm)','fontsize',14)
ylabel('t (ms)','fontsize',14)
zlabel('v (mV)','fontsize',14)

hold off


