

nlevels = 8; % number of grid levels. 1 means no multigrid, 2 means one coarse grid. etc
NX = 1*2^(nlevels-1); % Nx and Ny are given as function of grid levels
NY = 1*2^(nlevels-1);

uann = zeros(NX+2, NY+2); 
DX = 1.0/NX;
DY = 1.0/NY;
%disp(NX)

xc = linspace(0.5*DX, 1-0.5*DX, NX);
yc = linspace(0.5*DY, 1-0.5*DY, NY);
%disp(xc)
[XX, YY] = meshgrid(xc, yc);
exact = exact_solution(XX,YY);

surf(XX, YY,exact, 'EdgeColor', 'red')
xlabel('X');
ylabel('Y');
zlabel('U');
