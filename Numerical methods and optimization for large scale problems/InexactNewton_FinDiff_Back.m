function [Xk_IN, F_k_IN, G_k_norm_IN, k_IN, Xseq_IN, btseq_IN] = InexactNewton_FinDiff_Back...
    (X, k_max, tolgrad, c1, rho, bt_max, N, NEXT, h, FT, pcg_maxit)

% Description - Help of the function
%
% function [Xk, F_k, G_k_norm, k, Xseq, btseq] = InexactNewton_FinDiff_Back
%          (X, k_max, tolgrad, c1, rho, bt_max, N, NEXT, h, FT, pcg_maxit)        
%
% Function that performs Inexact Newton Method - Aspects to take into account:
% 1. Gradient and Hessian of F obtained by finite difference approximations
% 2. Preconditioned Conjugate Gradient Method to solve this: 
%    H_k * p + G_k <= eta_k * G_k
% 3. Newton Method complemented with a backtracking strategy (line search)
%
% INPUTS:
% X = n-dimensional initial column vector;
% k_max = maximum number of iterations permitted;
% tolgrad = value used as stopping criterion w.r.t. the norm of the
% gradient;
% c1 = ﻿the factor of the Armijo condition that must be a scalar in (0,1);
% rho = ﻿fixed factor, lesser than 1, used for reducing alpha0;
% bt_max = ﻿maximum number of steps for updating alpha during the
% N = dimension of the problem;
% NEXT = number of the problem selected;
% h = approximation step used for the finite difference computation of G;
% FT = Forcing Terms - Possible values 
%           FT = 1 - Linear rate convergence
%           FT = 2 - Superlinear rate convergence
%           FT = 3 - Quadratic rate convergence
% pcg_maxit = maximum number of iterations for the pcg solver
%
%
% OUTPUTS:
% Xk = last x computed by the function (best solution);
% F_k = value of F(Xk);
% G_k_norm = value of the norm of G(Xk)
% k = index of the last iteration performed
% Xseq = Matrix (NxK) with all solutions Xk computed during the iterations
% btseq = Row vector (size k) where elements are the number of backtracking
% iterations at each optimization step.


% Armijo Condition - Function handle
farmijo = @(F_k, alpha, G_k, pk) F_k + c1 * alpha * G_k' * pk;

% Initializations
Xseq_IN = zeros(N, k_max);
btseq_IN = zeros(1, k_max);

Xk_IN = X;
F_k_IN = TFFU28 (N, Xk_IN, NEXT);
k_IN = 0;
[G_k, H_k] = TFGHU28(N, Xk_IN, NEXT, h);
G_k_norm_IN = norm(G_k);
NonConvex = 0;

while k_IN < k_max && G_k_norm_IN >= tolgrad && NonConvex == 0
    
    switch FT
        case 1
            % Rate convergence is linear
            eta_k = 0.5;
            
        case 2
            % Rate convergence is superlinear
            eta_k = min(0.5, sqrt(G_k_norm_IN));
            
        case 3
            % Rate convergence is quadratic
            eta_k = min(0.5, G_k_norm_IN);
    end
            
    % INEXACT NEWTON METHOD
    % Not solving this linear system: H_k * p = - G_k
    % Instead of this, we approximate H_k * p + G_k <= small quantity
    % Small quantity: eta_k * G_k 
    
    % The smaller G_k is, the more precise is the direction p (because we 
    % are close to the solution)
    
    % Preconditioned Conjugate Gradient Method to solve this: 
    % H_k * p + G_k <= eta_k * G_k
    R = ichol(sparse(H_k));
    pk = pcg(H_k, -G_k, eta_k, pcg_maxit, R, R');
    
    % Reset alpha value
    alpha = 1;
    
    % Column vector - New candidate Xk
    Xnew = Xk_IN + alpha * pk;
    % Compute the value of f in the candidate new xk
    F_new = TFFU28(N,Xnew,NEXT);
    
    % Application of Backtracking strategy
    bt = 0;  
    while bt < bt_max && F_new > farmijo(F_k_IN, alpha, G_k, pk)
        % Reduce the value of alpha 
        alpha = rho * alpha;
        % Update X_new and F_new w.r.t. the reduced alpha
        Xnew = Xk_IN + alpha * pk;
        F_new = TFFU28(N,Xnew,NEXT);
        % Increase bactracking counter ("inner iteration")
        bt = bt + 1;
    end
    
    % Update Xk, F_k, G_k
    Xk_IN = Xnew;
    F_k_IN = F_new;
    [G_k, H_k] = TFGHU28(N,Xk_IN,NEXT,h);
    G_k_norm_IN = norm(G_k);
    
    % Increase step counter ("outer iteration")
    k_IN = k_IN + 1;
    
    % Store current xk in xseq
    Xseq_IN(:, k_IN) = Xk_IN;
    % Store bt iterations in btseq
    btseq_IN(k_IN) = bt;
    
    M = eig (H_k);
    
    for i = 1:N
        if M(i) < 0
            NonConvex = 1;
        end
    end
end

% "Cut" xseq and btseq to the correct size
Xseq_IN = [X Xseq_IN(:, 1:k_IN)];
btseq_IN = btseq_IN(1:k_IN);

end