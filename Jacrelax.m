function [u, res] = Jacrelax(level, nx, ny, u, f, iters, pre)
pre = false;
dx = 1.0/nx; dy = 1.0/ny;
Ax = 1.0/dx^2; Ay = 1.0/dy^2;
Ap = 1.0/(2.0*(Ax+Ay));

% Dirichlet BC
u(1,:) = -u(2,:);
u(end,:) = -u(end-1,:);
u(:,1) = -u(:,2);
u(:,end) = -u(:,end-1);

if pre && level > 1
    u(2:nx+1,2:ny+1) = -Ap*f(2:nx+1,2:ny+1);
    % Dirichlet BC
    u(1,:) = -u(2,:);
    u(end,:) = -u(end-1,:);
    u(:,1) = -u(:,2);
    u(:,end) = -u(:,end-1);

    iters = iters-1;
end

  for it = 0:iters+2
      for i = 2:nx+1
          for j = 2:ny+1
            u(i,j) = Ap*(Ax*(u(i+1,j) + u(i-1,j)) ...
            + Ay*(u(i,j+1) + u(i,j-1)) ...
            - f(i,j));
          end
      end
    % Dirichlet BC
        u(1,:) = -u(2,:);
        u(end,:) = -u(end-1,:);
        u(:,1) = -u(:,2);
        u(:,end) = -u(:,end-1);
  end
 res = zeros(nx+2, ny+2);
for i=2:nx+1
    for j=2:ny+1
        res(i,j) = f(i,j) - ((Ax*(u(i+1,j)+u(i-1,j)) + Ay*(u(i,j+1)+u(i,j-1)) - 2.0*(Ax+Ay)*u(i,j)));
    end
end


end
