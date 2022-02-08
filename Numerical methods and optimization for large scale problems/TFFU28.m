function [F] = TFFU28 (N, X, NEXT)

%
% Description - Help of the function
%
% function [F] = TFFU28(N, X, NEXT)
%
% Function that computes the value of different functions (from several
% selected problems) in a point X.
%
% INPUTS:
% N = dimension of the problem
% X = N-dimensional column vector;
% NEXT = number of the problem selected
%
% OUTPUTS:
% F = scalar value corresponding to the function evaluation in point X
% (depending on NEXT, the problem selected)

switch NEXT
    
    case 0
        F=(100*(X(2)-(X(1))^2)^2)+(1-X(1))^2;
        
    case 1
       aux = 0;
       for i=2:N
           suma = aux + (100*((X(i-1)^2-X(i))^2)+((X(i-1)-1)^2));
           aux = suma;
       end
       F = suma;
       
    case 2
       k = (N-2)/2;
       aux = 0;
       for j = 1:k
           i = 2*j;
           suma = aux + 100*((X(i-1)^2-X(i))^2)+(X(i-1)-1)^2+90*((X(i+1)^2-X(i+2))^2)+...
                      +(X(i+1)-1)^2+10*((X(i)+X(i+2)-2)^2)+((X(i)-X(i+2))^2)/10;
           aux = suma;
       end
       F = suma;
       
    case 3
       k = (N-2)/2;
       aux = 0;
       for j = 1:k
           i = 2*j;
           suma = aux + (X(i-1)+10*X(i))^2 + 5*((X(i+1)-X(i+2))^2)+...
                      +(X(i)-2*X(i+1))^4 + 10*(X(i-1)-X(i+2))^4;
           aux = suma;
       end
       F = suma;
end
       
       
        