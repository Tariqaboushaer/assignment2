clear all
clearvars
close all
clc
set(0, 'DefaultFigureWindowStyle', 'docked')

%Tariq Aboushaer
%101064544

%% QUESTION 2

% Solving current flow using Finite Diffrenece Method.

nx = 75;
ny = 50;
L = 20;
W = 10;
V1 = 1; 


SigmaC = 1;
SigmaI = 10e-2;
CM = SigmaC*ones(nx, ny);
CM(1:W,(1:L)+ny/2-L/2) = SigmaI;
CM((1:W)+nx-W,(1:L)+ny/2-L/2) = SigmaI;

figure();
hold on;
surf(linspace(0,1.5,ny), linspace(0,1,nx), CM,'EdgeColor','none','LineStyle','none');
xlabel('x');
ylabel('y');
zlabel('Conduction (Mho)');
view([120 25])
title('Sigma(x,y) (TA 101064544)')

%% PART (a)

V = NumSoln(nx, ny, CM, Inf, Inf, 0, V1);

figure();
hold on;
surf(linspace(0,1.5,ny), linspace(0,1,nx), V,'EdgeColor','none','LineStyle','none');
xlabel('x');
ylabel('y');
zlabel('Voltage (V)');
view([120 25])
colorbar
title('V(x,y) (TA 101064544)')

[Ex, Ey] = gradient(V);
Ex = -Ex;
Ey = -Ey;

figure();
quiver(linspace(0,1.5,ny), linspace(0,1,nx), Ex, Ey);
ylim([0 1]);
xlim([0 1.5]);
xlabel('x');
ylabel('y');
title('E(x,y) (TA 101064544)')

Jx = CM.*Ex;
Jy = CM.*Ey;
J = sqrt(Jx.^2 + Jy.^2);

figure();
hold on;
contourf(linspace(0,1.5,ny), linspace(0,1,nx), J,'EdgeColor','none','LineStyle','none');
quiver(linspace(0,1.5,ny), linspace(0,1,nx), Jx, Jy);
xlabel('x');
ylabel('y');
colorbar
title('J(x,y) (TA 101064544)')
%% PART (b)

figure();
hold on;
range = 20:5:100;
I = [];
for x = range
    
 I = [I TotalCurrent(x, ny, V1, SigmaC, SigmaI, W, L)];
    
end
plot(range, I);
ylabel('Current (A)');
xlabel('Mesh size');
title('Total current vs Width Mesh Size (TA 101064544)')
%% PART (c)

figure();
range = 0:1:50;
I = [];
for w = range
    
    I = [I TotalCurrent(nx, ny, V1, SigmaC, SigmaI, w, L)];
end
plot(range, I);
ylabel('Current (A)');
xlabel('Width');
title('Total current vs Box Width (TA 101064544)')

%% PART (d)

figure();
hold on;
range = logspace(-5,0, 50);
I = [];
for sigma = range
    
    I = [I TotalCurrent(nx, ny, V1, SigmaC, sigma, W, L)];

end
plot(range, I);
ylabel('Current (A)');
xlabel('Conduction (Mho)');
title('Total current vs Box Conduction (TA 101064544)')

%% Functions

function I = TotalCurrent(nx, ny, V1, SigmaC, SigmaI, W, L)
%total current
    CM = SigmaC*ones(nx, ny);
    CM(1:W,(1:L)+ny/2-L/2) = SigmaI;
    CM((1:W)+nx-W,(1:L)+ny/2-L/2) = SigmaI;
    V = NumSoln(nx, ny, CM, Inf, Inf, 0, V1);
    [Ex, Ey] = gradient(V);
    Ex = -Ex;
    Ey = -Ey;
    Jx = CM.*Ex;
    I = (abs(sum(Jx(1,:))) + abs(sum(Jx(nx,:))))/2;
end
