function [Xk_SD, F_k_SD, G_k_norm_SD, k_SD, Xseq_SD, btseq_SD] = ...
    SD_FinDiff_Back(X, k_max, tolgrad, c1, rho, bt_max, N, NEXT, h)
%
% [xk, fk, gradfk_norm, k, xseq] = steepest_descent(x0, f, gradf, alpha, kmax,
% tollgrad)
%
% Function that performs the steepest descent optimization method, for a 
% given function for the choice of the step length alpha.
%
% INPUTS:
% x0 = n-dimensional column vector;
% f = function handle that describes a function R^n->R;
% gradf = function handle that describes the gradient of f;
% alpha0 = the initial factor that multiplies the descent direction at each
% iteration;
% kmax = maximum number of iterations permitted;
% tolgrad = value used as stopping criterion w.r.t. the norm of the
% gradient;
% c1 = ﻿the factor of the Armijo condition that must be a scalar in (0,1);
% rho = ﻿fixed factor, lesser than 1, used for reducing alpha0;
% btmax = ﻿maximum number of steps for updating alpha during the 
% backtracking strategy.
%
% OUTPUTS:
% xk = the last x computed by the function;
% fk = the value f(xk);
% gradfk_norm = value of the norm of gradf(xk)
% k = index of the last iteration performed
% xseq = n-by-k matrix where the columns are the xk computed during the 
% iterations
% btseq = 1-by-k vector where elements are the number of backtracking
% iterations at each optimization step.
%

% Armijo Condition - Function handle
farmijo = @(F_k, alpha, G_k, pk) F_k + c1 * alpha * G_k' * pk;

% Initializations
Xseq_SD = zeros(N, k_max);
btseq_SD = zeros(1, k_max);

Xk_SD = X;
F_k_SD = TFFU28 (N, Xk_SD, NEXT);
k_SD = 0;
[G_k, ~] = TFGHU28(N, Xk_SD, NEXT, h);
G_k_norm_SD = norm(G_k);
NonConvex = 0;
alpha0 = 5;

while k_SD < k_max && G_k_norm_SD >= tolgrad && NonConvex == 0
    % Compute the descent direction
    pk = -G_k;
    
    % Reset the value of alpha
    alpha = alpha0;
    
    % Compute the candidate new xk
    Xk_new = Xk_SD + alpha * pk;
    % Compute the value of f in the candidate new xk
    Fk_new = TFFU28 (N, Xk_new, NEXT);
    
    bt = 0;
    % Backtracking strategy: 
    % 2nd condition is the Armijo condition not satisfied
    while bt < bt_max && Fk_new > farmijo(F_k_SD, alpha, G_k, pk)
        % Reduce the value of alpha
        alpha = rho * alpha;
        % Update xnew and fnew w.r.t. the reduced alpha
        Xk_new = Xk_SD + alpha * pk;
        Fk_new = TFFU28 (N, Xk_new, NEXT);
        
        % Increase the counter by one
        bt = bt + 1;
        
    end
    
    % Update xk, fk, gradfk_norm
    Xk_SD = Xk_new;
    F_k_SD = Fk_new;
    [G_k, ~] = TFGHU28(N, Xk_SD, NEXT, h);
    G_k_norm_SD = norm(G_k);
    
    % Increase the step by one
    k_SD = k_SD + 1;
    
    % Store current xk in xseq
    Xseq_SD(:, k_SD) = Xk_SD;
    % Store bt iterations in btseq
    btseq_SD(k_SD) = bt;
end

% "Cut" xseq and btseq to the correct size
Xseq_SD = Xseq_SD(:, 1:k_SD);
btseq_SD = btseq_SD(1:k_SD);

end