function [X, IERR, FMIN, XMAX] = TIUD28 (N, NEXT)

%
% Description - Help of the function
%
% function [X, IERR, FMIN, XMAX] = TIUD28(N, NEXT)
%
% Function that computes the initial X in different problems.
%
% INPUTS:
% N = dimension of the problem
% NEXT = number of the problem selected
%
% OUTPUTS:
% X = column N-dimensional vector 
% IERR = error indicator (0 - correct data; 1 - N is too small)
% FMIN = lower bound of the objective function value
% XMAX = maximum stepsize

X = zeros(N,1);

switch NEXT
 
    case 1
        FMIN = 0;
        XMAX = 1;
        if N <= 1
            IERR = 1;
        else
            for i=1:N
                if mod(i,2)==1
                    X(i) = -1.2;
                else
                    X(i) = 1;
                end
            end
            IERR = 0;
        end
        
        
    case 2
        FMIN = 0;
        XMAX = 1;
        if N <= 1
            IERR = 1;
        else
            for i = 1:N
                if i<=4
                    if mod(i,2)==1
                        X(i) = -3;
                    else
                        X(i) = -1;
                    end
                else
                    if mod(i,2)==1
                        X(i) = -2;
                    else
                        X(i) = 0;
                    end
                end
            end
            IERR = 0;
        end
        
        
    case 3
        FMIN = 0;
        XMAX = 1;
        if N <= 1
            IERR = 1;
        else
            for i = 1:N
                if mod(i,4)==1
                    X(i) = 3;
                else if mod(i,4)==2
                    X(i) = -1;
                else if mod(i,4)==3
                    X(i) = 0;
                else
                    X(i) = 1;
                end
                end
                end
            end
            IERR = 0;
        end                
end



