%% EVALUATION CODE - ROSENBROCK FUNCTION - N = 2

N = 2;
NEXT = 0;
X1 = [1.2; 1.2];
X2 = [-1.2; 1];
k_max = 1e4;
tolgrad = 1e-6;
c1 = 1e-4;
rho = 0.5;
bt_max = 100;
h = 1e-8;
FT = 3;
pcg_maxit = 100*N;

% Nelder Method
k_max_ND = 70;
rho_ND = 1;
chi=2;
gamma=0.5;
sigma=0.5;

%% Newton Method 
tic
[Xk_N1, F_k_N1, G_k_norm_N1, k_N1, Xseq_N1, btseq_N1] = ...
    Newton_FinDiff_Back(X1, k_max, tolgrad, c1, rho, bt_max, N, NEXT, h);
toc

tic
[Xk_N2, F_k_N2, G_k_norm_N2, k_N2, Xseq_N2, btseq_N2] = ...
    Newton_FinDiff_Back(X2, k_max, tolgrad, c1, rho, bt_max, N, NEXT, h);
toc

%% Inexact Newton Method
tic
[Xk_IN1, F_k_IN1, G_k_norm_IN1, k_IN1, Xseq_IN1, btseq_IN1] = InexactNewton_FinDiff_Back...
    (X1, k_max, tolgrad, c1, rho, bt_max, N, NEXT, h, FT, pcg_maxit);
toc

tic
[Xk_IN2, F_k_IN2, G_k_norm_IN2, k_IN2, Xseq_IN2, btseq_IN2] = InexactNewton_FinDiff_Back...
    (X2, k_max, tolgrad, c1, rho, bt_max, N, NEXT, h, FT, pcg_maxit);
toc

%% Steepest Descent Method
tic
[Xk_SD1, F_k_SD1, G_k_norm_SD1, k_SD1, Xseq_SD1, btseq_SD1] = ...
    SD_FinDiff_Back(X1, k_max, tolgrad, c1, rho, bt_max, N, NEXT, h);
toc

tic
[Xk_SD2, F_k_SD2, G_k_norm_SD2, k_SD2, Xseq_SD2, btseq_SD2] = ...
    SD_FinDiff_Back(X2, k_max, tolgrad, c1, rho, bt_max, N, NEXT, h);
toc

%% Nelder-Mead Method
tic
[X0_ND,Xk_ND, F_k_ND, k_ND, Xseq_ND] = Nelder_Method...
    (k_max_ND, N, NEXT,rho_ND,chi,gamma,sigma);
toc

%% PLOTS 

% Creation of the meshgrid for the contour-plot
[L, M] = meshgrid(linspace(-6, 6, 500), linspace(-6, 6, 500));
% Computation of the values of f for each point of the mesh
Z = 100*(M-L.^2).^2+(1-L).^2;

% PLOT - NEWTON METHOD WITH FINITE DIFFERENCES

fig1 = figure();
sgtitle ('NEWTON METHOD WITH FD - - TWO STARTING POINTS')
subplot (2,1,1)
contour(L, M, Z);
hold on
title ('Contour Graph') , xlabel ('Eje x'), ylabel ('Eje y')
plot(Xseq_N1(1,:), Xseq_N1(2,:), '--om')
plot(Xseq_N2(1,:), Xseq_N2(2,:), '--ok')
plot(1,1,'*r')
hold off

subplot (2,1,2)
meshc(L, M, Z);
hold on
title ('Mesh and contour graph') , xlabel ('Eje x'), ylabel ('Eje y'), zlabel ('Eje z')
plot(Xseq_N1(1,:), Xseq_N1(2,:), '--om')
plot(Xseq_N2(1,:), Xseq_N2(2,:), '--ok')
plot(1,1,'*r')
hold off

% PLOT - INEXACT NEWTON METHOD WITH FINITE DIFFERENCES

fig2 = figure();
sgtitle ('INEXACT NEWTON METHOD WITH FD - TWO STARTING POINTS')
subplot (2,1,1)
contour(L, M, Z);
hold on
title ('Contour Graph') , xlabel ('Eje x'), ylabel ('Eje y')
plot(Xseq_IN1(1,:), Xseq_IN1(2,:), '--om')
plot(Xseq_IN2(1,:), Xseq_IN2(2,:), '--ok')
plot(1,1,'*r')
hold off

subplot (2,1,2)
meshc(L, M, Z);
hold on
title ('Mesh and contour graph') , xlabel ('Eje x'), ylabel ('Eje y'), zlabel ('Eje z')
plot(Xseq_IN1(1,:), Xseq_IN1(2,:), '--om')
plot(Xseq_IN2(1,:), Xseq_IN2(2,:), '--ok')
plot(1,1,'*r')
hold off

% PLOT - STEEPEST DESCENT WITH FINITE DIFFERENCES

fig3 = figure();
sgtitle ('STEEPEST METHOD WITH FD - TWO STARTING POINTS')
subplot (2,1,1)
contour(L, M, Z);
hold on
title ('Contour Graph') , xlabel ('Eje x'), ylabel ('Eje y')
plot(Xseq_SD1(1,:), Xseq_SD2(2,:), '--om')
plot(Xseq_SD2(1,:), Xseq_SD2(2,:), '--ok')
plot(1,1,'*r')
hold off

subplot (2,1,2)
meshc(L, M, Z);
hold on
title ('Mesh and contour graph') , xlabel ('Eje x'), ylabel ('Eje y'), zlabel ('Eje z')
plot(Xseq_SD1(1,:), Xseq_SD1(2,:), '--om')
plot(Xseq_SD2(1,:), Xseq_SD2(2,:), '--ok')
plot(1,1,'*r')
hold off

%PLOT - NELDER-MEAD METHOD

fig4 = figure();
sgtitle ('NELDER-MEAD METHOD ')
subplot (2,1,1)
contour(L, M, Z);
hold on
title ('Contour Graph (First Iteration)') , xlabel ('Eje x'), ylabel ('Eje y')

plot(Xseq_ND(1,1:4), Xseq_ND(2,1:4), '--ob')
plot(Xseq_ND(1,29:32), Xseq_ND(2,29:32), '--ok')
plot(Xseq_ND(1,41:44), Xseq_ND(2,41:44), '--og')

b1=plot(Xk_ND(1),Xk_ND(2),'hk');
b2=plot(1,1,'*r');
legend([b1 b2],'Nelder-Mead Optimal','Real Exact Optimal')

hold off

subplot (2,1,2)
contour(L, M, Z);
hold on
title ('Contour Graph Zoomed In  (Last Iteration)') , xlabel ('Eje x'), ylabel ('Eje y')

plot(Xseq_ND(1,1:4), Xseq_ND(2,1:4), '--ob')
plot(Xseq_ND(1,29:32), Xseq_ND(2,29:32), '--ok')
plot(Xseq_ND(1,41:44), Xseq_ND(2,41:44), '--og')

b1=plot(Xk_ND(1),Xk_ND(2),'hk');
b2=plot(1,1,'*r');
legend([b1 b2],'Nelder-Mead Optimal','Real Exact Optimal')

hold off
