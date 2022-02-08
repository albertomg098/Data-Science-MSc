function [G, H] = TFGHU28(N,X,NEXT,h)

% Description - Help of the function
%
% function [G,H] = TFGHU28(N, X, NEXT)
%
% Function that approximate the Gradient and the Hessian of F(evaluated in 
% TFFU28) in X (column vector) with the forward finite difference method.
%
% INPUTS:
% N = dimension of the problem
% X = N-dimensional column vector;
% h = approximation step used for the finite difference computation of G
% NEXT = number of the problem selected
%
% OUTPUTS:
% G = column vector (same size of X) corresponding to the approximation
% of the Gradient of F in X.
% H = matrix (NxN) corresponding to the approximation of the Hessian of F
% in X.

G = zeros (N,1);
        
    for i=1:N
        X_h = X;
        X_h(i) = X_h(i) + h;
        G(i) = (TFFU28 (N, X_h, NEXT) - TFFU28 (N, X, NEXT))/ h;
    end
    
H = zeros(N, N);
% h for the Hessian is suggested to be greater the the h used for gradient
% in order to avoid numericall cancelation problems (h_hess = sqrt(h_grad))
h1 = sqrt(h);

for j=1:N
    % Elements on the diagonal
    Xh_plus = X;
    Xh_minus = X;
    Xh_plus(j) = Xh_plus(j) + h1;
    Xh_minus(j) = Xh_minus(j) - h1;
    H(j,j) = (TFFU28 (N, Xh_plus, NEXT) - 2*TFFU28 (N, X, NEXT) + TFFU28 (N, Xh_minus, NEXT))/(h1^2);
    
    % Elements of the other elements
    for i=j+1:N
        Xh_plus_ij = X;
        Xh_plus_ij([i, j]) = Xh_plus_ij([i, j]) + h1;
        Xh_plus_i = X;
        Xh_plus_i(i) = Xh_plus_i(i) + h1;
        Xh_plus_j = X;
        Xh_plus_j(j) = Xh_plus_j(j) + h1;
        H(i,j) = (TFFU28(N,Xh_plus_ij,NEXT) - TFFU28(N,Xh_plus_i,NEXT) - TFFU28(N,Xh_plus_j,NEXT) +...
                  TFFU28(N,X,NEXT))/(h1^2);
        H(j,i)=H(i,j);
    end
end

end