%% EVALUATION OF THE CODE - N-DIMENSIONAL PROBLEMS

% Initializions
N = 500;
NEXT = 2;
[X, IERR, FMIN, XMAX] = TIUD28 (N, NEXT);
k_max = 2*N;
tolgrad = 1e-4;
c1 = 1e-4;
rho = 0.5;
bt_max = 100;
h = 1e-8;

% Nelder Method
k_max_ND = 1000;
rho_ND = 1;
chi=2;
gamma=0.5;
sigma=0.5;

%% Steepest Descent Method
tic
[Xk_SD, F_k_SD, G_k_norm_SD, k_SD, Xseq_SD, btseq_SD] = ...
    SD_FinDiff_Back(X, k_max, tolgrad, c1, rho, bt_max, N, NEXT, h);
toc

%% NELDER-MEAD method
tic
[X0,Xk_ND, F_k_ND, k_ND] = Nelder_Method(k_max_ND, N, NEXT, rho_ND, chi, gamma, sigma);
toc

