function [X0,Xk_N, F_k_N, k_N,Xseq_N] = ...
    Nelder_Method(k_max, N, NEXT,rho,chi,gamma,sigma)

% Description - Help of the function
%
% function [X0,Xk_N, F_k_N, k_N,Xseq_N] = ...
%           Nelder_Method(k_max, N, NEXT,rho,chi,gamma,sigma)
%
% INPUTS:
% k_max = maximum number of iterations permitted;
% N = dimension of the problem;
% NEXT = number of the problem selected;
% tolgrad = value used as stopping criterion w.r.t. the norm of the
% gradient;
% c1 = ﻿the factor of the Armijo condition that must be a scalar in (0,1);
%
% Typical values:
%  rho=1
%  chi= 2
%  gamma=0.5
%  sigma=0.5
%
% OUTPUTS:
% X0 = initial simplex set (randomly selected in N-Dimension);
% Xk_N = final optimal solution obtained by Nelder Mead Method;
% F_k_N = final optimum (value of the function at the optimal point);
% Xseq_N = When NEXT = 0, it represents the group of simplex sets;
%
%
% Initialize the matrix for defining the simplex(each column is a vertice of
% the simple)

switch NEXT
    case 0
        X=[1.2,-1.2,0;1.2,1,0];
        %Generate randomly a third point with values between 0-2
        for i=1:2
            X(i,3)=2*rand();
        end
        Xseq_N=[X X(:,1)];
        
    otherwise
        X=zeros(N,N+1);
        %Generate randomly points with values between 0-2
        for i=1:N
            for j=1:(N+1)
                X(i,j)=2*rand();
            end
        end
end

X0=X;
F=zeros(1,N+1);
k=0;

%Start the loop

while k<k_max 
    k=k+1;
    
%1. ORDERING PHASE
%Compute the value of the function at he vertices
for i=1:N+1
    F(i)=TFFU28(N,X(:,i),NEXT);
end

%Order both the function vector and the points
[F_ordered,I]=sort(F);

X_ordered=zeros(N,N+1);
for j=1:(N+1)
    X_ordered(:,j)=X(:,I(j));
end

%2. REFLECTION PHASE
%Compute the barycenter

x_baricenter=zeros(N,1);
for j=1:N
    x_baricenter=x_baricenter+X_ordered(:,j);
end
x_baricenter=x_baricenter/N;

%Compute the reflection of x(N+1)
X_R=zeros(N,1);
X_E=zeros(N,1);
X_C=zeros(N,1);

X_R=x_baricenter+(rho*(x_baricenter-X_ordered(:,N+1)));
F_R=TFFU28(N,X_R,NEXT);

if F_R>=F_ordered(1) && F_R<F_ordered(N)
    %Enought reduction with xR%
    %Accept xR and go to next step
    X_ordered(:,N+1)=X_R;
else
    if F_R<F_ordered(1)
        %Large reduction 
        X_E=x_baricenter+(chi*(X_R-x_baricenter));
        F_E=TFFU28(N,X_E,NEXT);
        %EXPANSION PHASE
        if F_E<F_R
            %Accept xE and go to next step
            X_ordered(:,N+1)=X_E;
        else
            %Accept xR and go to next step
            X_ordered(:,N+1)=X_R;
        end
    else
        %Poor Reduction (fR>=fN)
        %CONTRACTION PHASE%
        X_C=x_baricenter+(gamma*(X_ordered(:,N+1)-x_baricenter));
        F_C= TFFU28(N,X_C,NEXT);
        if F_C<F_ordered(N+1)
            %Accept xC and go to next step
            X_ordered(:,N+1)=X_C;
        else
            %SHRINKAGE PHASE
            for i=2:(N+1)
                X_ordered(:,i)=X_ordered(:,1)-(sigma*(X_ordered(:,i)-X_ordered(:,1)));
            end
        end
    end
end

%Determine Xseq for base case(N=2)
if NEXT==0
    g=k_max/10;
    if mod(k,g)==0
        Xseq_N=[Xseq_N X_ordered X_ordered(:,1)];
    end
end

%Reset values for the next iteration
X=X_ordered;

end

k_N=k;
%Order the solutions to report the best one
for i=1:N+1
    F(i)=TFFU28(N,X(:,i),NEXT);
end

[F_ordered,I]=sort(F);

X_ordered=zeros(N,N+1);
for j=1:(N+1)
    X_ordered(:,j)=X(:,I(j));
end

Xk_N=X_ordered(:,1);
F_k_N=F_ordered(1);
end