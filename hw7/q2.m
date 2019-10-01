% MATH 151B, HOMEWORK 7, Question 2
% WANG, ZHENG (404855295)
A = [3 3 3; 4 9 2; 5 2 3];
q_1 = 5;
q_2 = 2;
x = [1;1;1];

% q = 5
fprintf('Using q = 5\n')
[lambda, v] = power_method(A,x,q_1,10^-5,10000);
fprintf('The eigenvalue find is:\n')
disp(lambda)
fprintf('The corresponding eigenvector find is:\n')
disp(v)

% q = 2
fprintf('Using q = 2\n')
[lambda, v] = power_method(A,x,q_2,10^-5,10000);
fprintf('The eigenvalue find is:\n')
disp(lambda)
fprintf('The corresponding eigenvector find is:\n')
disp(v)

% INPUTS:
% A - the matrix to be solved
% q - Find eigenvalue closest to q
% x - The initial vector
% Tol - Tolerance
% N - max iteration number
function [lambda, v] = power_method(A,x,q,Tol,N)
    n = size(A,1);
    M = A - q*eye(n);
    k = 1;
    [~, p] = max(abs(x));  
    while k<=N
        y = M\x;
        % mu is the estimate of eigenvalue
        % p is the position of the largest entry
        mu = y(p);
        % update p
        [~, p] = max(abs(y));  
        yp = y(p);
        err = max(abs(x-y/yp));
        x = y/yp;
        if err < Tol
            lambda = 1/mu + q;
            v = x;
            return;
        end
        k = k+1;
    end
    lambda = mu;
    v = x;
    fprintf('Reach max iteration')
end
        
    