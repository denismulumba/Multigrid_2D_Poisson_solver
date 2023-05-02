% Set up the problem parameters

L = 1;      % length of the domain in x and y directions
Nx = 128;    % number of grid points in x direction
Ny = 128;    % number of grid points in y direction
dx = L/(Nx);  % grid spacing in x direction
dy = L/(Ny);  % grid spacing in y direction

% Define the boundary conditions
u0 = 0;     % Dirichlet boundary condition on u(x,y) at x=0
uL = 0;     % Dirichlet boundary condition on u(x,y) at x=L
v0 = 0;     % Dirichlet boundary condition on u(x,y) at y=0
vL = 0;     % Dirichlet boundary condition on u(x,y) at y=L

% Define the forcing function
f = @(x,y) -2*pi^2*sin(pi*x).*sin(pi*y);  % Laplace operator of sin(pi*x)*sin(pi*y)

% Set up the grid
x = linspace(dx,L-dx,Nx);
y = linspace(dy,L-dy,Ny);
[X,Y] = meshgrid(x,y);

% Initialize the solution
U = zeros(Nx,Ny);

% Set up the finite difference stencil
a = 1/dx^2;
b = 1/dy^2;
c = -2*(a+b);

% Set up the coefficient matrix
tb = tic;
A = zeros(Nx*Ny,Nx*Ny);
for i = 1:Nx
    for j = 1:Ny
        k = (i-1)*Ny + j;
        if i == 1
            A(k,k) = 1;
        elseif i == Nx
            A(k,k) = 1;
        elseif j == 1
            A(k,k) = 1;
        elseif j == Ny
            A(k,k) = 1;
        else
            A(k,k-Ny) = a;
            A(k,k-1) = b;
            A(k,k) = c;
            A(k,k+1) = b;
            A(k,k+Ny) = a;
        end
    end
end

% Set up the right-hand side
b = zeros(Nx*Ny,1);
for i = 1:Nx
    for j = 1:Ny
        k = (i-1)*Ny + j;
        if i == 1
            b(k) = u0;
        elseif i == Nx
            b(k) = uL;
        elseif j == 1
            b(k) = v0;
        elseif j == Ny
            b(k) = vL;
        else
            b(k) = dx^2*f(x(i),y(j));
        end
    end
end

% Solve the linear system using backslash operator
Uvec = A\b;
y = toc(tb);
disp(y)
U = reshape(Uvec,Nx,Ny);

% Plot the solution
surf(X,Y,U);
xlabel('x');
ylabel('y');
zlabel('u(x,y)');
title('Solution to 2D Poisson Equation');
