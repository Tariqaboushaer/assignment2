clear all
clearvars
close all
clc
set(0, 'DefaultFigureWindowStyle', 'docked')

%Tariq Aboushaer
%101064544

%% QUESTION 1
% 2D Laplace case

W = 2;
L = 3;
V0 = 1;

dx = 0.2;
dy = 0.2;
nx = L/dx;
ny = W/dy;



C1 = -2*(1/dx^2 + 1/dy^2);
C2 = 1/(dx^2);
C3 = 1/(dy^2);
g = zeros(nx*ny,nx*ny);
f = zeros(nx*ny,1);

%% 
for x=2:(nx-1)
    for y=2:(ny-1)
        c = (y-1).*nx + x;
        g(c,c) = C1;
        g(c,((y-1).*nx + x-1)) = C2;
        g(c,((y-1).*nx + x+1)) = C2;
        g(c,((y-2).*nx + x)) = C3;
        g(c,((y).*nx + x)) = C3;
    end
end


for y=1:ny
    c = ((y-1).*nx + 1);
    g(c,c) = 1;
    
    f(c) = V0;
    
    c = ((y-1).*nx + nx);
    g(c,c) = 1;
end


for x=2:(nx-1)
    c = (1-1).*nx + x;
    g(c,c) = 1;
    g(c,(2-1).*nx + x) = -1;
    
    c = ((ny-1).*nx + x);
    g(c,c) = 1;
    g(c,((ny-2).*nx + x)) = -1;
end
%%
v = g\f;
v = reshape(v,[],ny)';

figure();
surf(linspace(0,L,nx),linspace(0,W,ny),v);
xlabel('x');
ylabel('y');
title(sprintf('2-D plot of V(x) = %.2f (TA 101064544)', dx));
set(gca, 'View', [45 45])

figure();
surf(linspace(0,L,nx),linspace(0,W,ny),v);
xlabel('x');
ylabel('y');
title(sprintf('2-D plot of V(x) = %.2f (TA 101064544)', dx));
set(gca, 'View', [45 45])
shading interp


Soln = zeros(ny, nx);
x1 = repmat(linspace(-L/2,L/2,nx),ny,1);
y1 = repmat(linspace(0,W,ny),nx,1)';
it = 100;
avgError = zeros(it,1);


for c=1:it
    n = 2*c - 1;
    Soln = Soln + 1./n.*cosh(n.*pi.*x1./W) ...
        ./cosh(n.*pi.*(L./2)./W).*sin(n.*pi.*y1./W);

    avgError(c) = mean(mean(abs(Soln.*4.*V0./pi - v)));
end

Soln = Soln.*4.*V0./pi;

figure();
surf(linspace(0,L,nx),linspace(0,W,ny),Soln);
shading interp
xlabel('x');
ylabel('y');
title(sprintf('Analytical Graph with %d iterations (TA 101064544)', it));


