% This is an example showing how to call the mgd2d solver.
%pre = false;
% analytical solution
Uann = @(x, y) sin(pi*x).*sin(pi*y);
%Uann = @(x,y) exp(x-y).*(x.^2 -1)*(y.^2 -1);
%Uann = @(x,y) (1/(5*pi^2))*sin(pi*x).*sin(2*pi*y);
%Uann = @(x,y) exp(sin(2*pi*(x+2*y)));

% RHS corresponding to above
source = @(x, y) -2*pi^2*sin(pi*x).*sin(pi*y);
%source = @(x,y) 4*(x.^2 + y.^2 - 1).*exp(x-y);
%source = @(x,y) sin(pi*x).*sin(2*pi*y);
%source = @(x,y) 10*pi^2*(1+cos(4*pi*(x+2*y))-2*sin(2*pi*(x+2*y))).*...
  %  exp(sin(2*pi*(x+2*y)));

% input
max_cycles = 3000; % maximum number of V cycles
nlevels = 9; % number of grid levels. 1 means no multigrid, 2 means one coarse grid. etc
NX = 1*2^(nlevels-1); % Nx and Ny are given as function of grid levels
NY = 1*2^(nlevels-1);
tol = 1e-9;
%omega = 4/5;

% the grid has one layer of ghost cells to help apply the boundary conditions
uann = zeros(NX+2, NY+2); % analytical solution
u = zeros(NX+2, NY+2); % approximation
f = zeros(NX+2, NY+2); % RHS
%disp(size(u))
% calculate the RHS and exact solution
DX = 1.0/NX;
DY = 1.0/NY;

disp(DX)
xc = linspace(0.5*DX, 1-0.5*DX, NX);
yc = linspace(0.5*DY, 1-0.5*DY, NY);
[XX, YY] = meshgrid(xc, yc);

uann(2:NX+1, 2:NY+1) = Uann(XX, YY);
f(2:NX+1, 2:NY+1) = source(XX, YY);

disp(['NX: ', num2str(NX), ', NY: ', num2str(NY), ', tol: ', num2str(tol), ', levels: ', num2str(nlevels)]);

% start solving
tb = tic;

% V cycle
for it = 1:max_cycles
    [u, res] = V_cycle(NX, NY, nlevels, u, f);
    rtol = max(max(abs(res)));
    if(rtol < tol)
            break;
    end
    error = uann(2:NX+1, 2:NY+1) - u(2:NX+1, 2:NY+1);
    disp(['  cycle: ', num2str(it), ', L_inf(res.) = ', num2str(rtol), ', L_inf(true error): ', num2str(max(max(abs(error))))]);
end

disp(['Elapsed time: ', num2str(toc(tb)), ' seconds']);

error = uann(2:NX+1, 2:NY+1) - u(2:NX+1, 2:NY+1);

disp(['L_inf (true error): ', num2str(max(max(abs(error))))]);

figure(1);
surf(XX, YY, abs(error))
%surf(XX, YY, u(2:NX+1, 2:NY+1), 'EdgeColor', 'blue');
xlabel('X');
ylabel('Y');
zlabel('U');
title('Approximation of U(x,y)');
