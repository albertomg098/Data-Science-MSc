function [Xk_N, F_k_N, G_k_norm_N, k_N, Xseq_N, btseq_N] = ...
    Newton_FinDiff_Back(X, k_max, tolgrad, c1, rho, bt_max, N, NEXT, h)

% Description - Help of the function
%
% function [Xk, F_k, G_k_norm, k, Xseq, btseq] = 
%   Newton_FinDiff_Back(X, k_max, tolgrad, c1, rho, bt_max, N, NEXT, h)        
%
% Function that performs Newton Method - Aspects to take into account:
% 1. Gradient and Hessian of F obtained by finite difference approximations 
% 2. Newton Method complemented with a backtracking strategy (line search)
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
Xseq_N = zeros(N, k_max);
btseq_N = zeros(1, k_max);

Xk_N = X;
F_k_N = TFFU28 (N, Xk_N, NEXT);
k_N = 0;
[G_k, H_k] = TFGHU28(N, Xk_N, NEXT, h);
G_k_norm_N = norm(G_k);
NonConvex = 0;


while k_N < k_max && G_k_norm_N >= tolgrad && NonConvex == 0
    
    % Preconditioned Conjugate Gradient Method to solve this: 
    % H_k * p + G_k <= eta_k * G_k
    pk = -H_k\G_k;

    % Reset alpha value
    alpha = 1;
    
    % Column vector - New candidate Xk
    Xnew = Xk_N + alpha * pk;
    % Compute the value of f in the candidate new xk
    F_new = TFFU28(N,Xnew,NEXT);
    
    % Application of Backtracking strategy
    bt = 0; 
    % While loop until Armijo condition is satisfied
    while bt < bt_max && F_new > farmijo(F_k_N, alpha, G_k, pk)
        % Reduce the value of alpha 
        alpha = rho * alpha;
        % Update X_new and F_new w.r.t. the reduced alpha
        Xnew = Xk_N + alpha * pk;
        F_new = TFFU28(N,Xnew,NEXT);
        % Increase bactracking counter ("inner iteration")
        bt = bt + 1;
    end
    
    % Update Xk, F_k, G_k
    Xk_N = Xnew;
    F_k_N = F_new;
    [G_k, H_k] = TFGHU28(N,Xk_N,NEXT,h);
    G_k_norm_N = norm(G_k);
    M = eig (H_k);
    
    for i=1:N
        if M(i)<=0
            NonConvex = 1;
        end
    end
   
    % Increase step counter ("outer iteration")
    k_N = k_N + 1;
    
    % Store current xk in xseq
    Xseq_N(:, k_N) = Xk_N;
    % Store bt iterations in btseq
    btseq_N(k_N) = bt;
end

% "Cut" xseq and btseq to the correct size
Xseq_N = [X Xseq_N(:, 1:k_N)];
btseq_N = btseq_N(1:k_N);

end