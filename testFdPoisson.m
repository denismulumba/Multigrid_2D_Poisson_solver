

function [u,x,y] = testFdPoisson(ffun, gfun, a, b, m)

h = (b-a)/(m+1);
[x,y] = meshgrid(a:h:b);

idx = 2:m+1;
idy = 2:m+1;

ubs = feval(gfun, x(1, 1:m+2), y(1, 1:m+2));

ubn = feval(gfun, x(m+2, 1:m+2), y(m+2, 1:m+2));

ube = feval(gfun, x(idy, m+2), y(idy, m+2));

ubw = feval(gfun, x(idx, 1), y(idx, 1));

f = feval(ffun, x(idy, idx), y(idy, idx));

f(:,1) = f(:,1) - ubw/h^2;
f(:,m) = f(:,m) - ube/h^2;
f(1,1:m) = f(1,1:m) - ubs(idx)/h^2;
f(m,1:m) = f(m, 1:m) - ubn(idx)/h^2;

f = reshape(f, m*m, 1);

D2 = toeplitz(spdiags([-2 1], [0,1], 1, m)); 
D2x = 1/h^2*kron(D2, speye(m));
D2y = 1/h^2*kron(speye(m), D2);

u = (D2x + D2y)\f;


u = reshape(u,m,m);

u = [ubs; [ubw, u, ube]; ubn];

end