


a = 0;
b = 1;
k = 8;
m = 2^k ;
h = (b-a)/(m+1);
%uexact = @(x,y) exp(sin(2*pi*(x+2*y)));
uexact = @(x, y) sin(pi*x).*sin(pi*y);

%f = @(x,y) 10*pi^2*(1+cos(4*pi*(x+2*y))-2*sin(2*pi*(x+2*y))).*...
    exp(sin(2*pi*(x+2*y)));
 f = @(x,y)  -2*pi^2*sin(pi*x).*sin(pi*y);

g = @(x,y) uexact(x,y);

tic
[u,x,y] = testFdPoisson(f, g, a, b, m);
gedirect = toc;

fprintf('Direct Gaussian elimination takes %d s\n', gedirect);

figure,
set(gcf, 'DefaultAxesFontSize',10,'PaperPosition', [0 0 3.5 3.5]),
surf(x, y, u),
xlabel('x')
ylabel('y')
zlabel('u(x,y)')
title(strcat('Numerical Solution to Poisson Equation, h =', num2str(h)));

%plot error
figure,
set(gcf, 'DefaultAxesFontSize',10,'PaperPosition', [0 0 3.5 3.5]),
surf(x, y, u-uexact(x,y)),
xlabel('x')
ylabel('y')
zlabel('Error')
title(strcat('Error, h =', num2str(h)));