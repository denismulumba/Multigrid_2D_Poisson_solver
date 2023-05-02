function [u, res] = Jacrelax(level, nx, ny, u, f, iters, pre)

% Set pre-processing flag to false
pre = false;

% Calculate grid spacings in x and y
dx = 1.0/nx;
dy = 1.0/ny;

% Calculate constants for discrete Laplacian operator
Ax = 1.0/dx^2;
Ay = 1.0/dy^2;
Ap = 1.0/(2.0*(Ax+Ay));

% Apply Dirichlet boundary conditions
u(1,:) = -u(2,:);
u(end,:) = -u(end-1,:);
u(:,1) = -u(:,2);
u(:,end) = -u(:,end-1);

% If pre-processing is enabled and this is not the coarsest level, 
% set u to the negative of the right-hand side
if pre && level > 1
    u(2:nx+1,2:ny+1) = -Ap*f(2:nx+1,2:ny+1);
    
    % Apply Dirichlet boundary conditions
    u(1,:) = -u(2,:);
    u(end,:) = -u(end-1,:);
    u(:,1) = -u(:,2);
    u(:,end) = -u(:,end-1);

    % Reduce the number of iterations to account for pre-processing
    iters = iters-1;
end

% Perform Jacobi relaxation for the specified number of iterations
for it = 0:iters+2
    for i = 2:nx+1
        for j = 2:ny+1
            % Update u using the discrete Laplacian operator and the right-hand side f
            u(i,j) = Ap*(Ax*(u(i+1,j) + u(i-1,j)) ...
            + Ay*(u(i,j+1) + u(i,j-1)) ...
            - f(i,j));
        end
    end
    
    % Apply Dirichlet boundary conditions
    u(1,:) = -u(2,:);
    u(end,:) = -u(end-1,:);
    u(:,1) = -u(:,2);
    u(:,end) = -u(:,end-1);
end

% Calculate the residual vector
res = zeros(nx+2, ny+2);
for i=2:nx+1
    for j=2:ny+1
        res(i,j) = f(i,j) - ((Ax*(u(i+1,j)+u(i-1,j)) + Ay*(u(i,j+1)+u(i,j-1)) - 2.0*(Ax+Ay)*u(i,j)));
    end
end

end
