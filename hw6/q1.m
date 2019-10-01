% MATH 151b, HW 6
% Question 1
fprintf('Solution with 4 iteration is:\n')
disp(Newton(@J,@F,4,[1;1;1]))
fprintf('Solution with 8 iteration is:\n')
disp(Newton(@J,@F,8,[1;1;1]))

% The function to solve
function Y = F(x1,x2,x3)
    y = zeros(3,1);
    y(1) = x1^2 + x2 -37;
    y(2) = x1 - x2^2 - 5;
    y(3) = x1 + x2 + x3 -3;
    Y = y;
end

% The Jacobian of the function
function Jac = J(x1,x2,x3)
    Jacobian = zeros(3);
    Jacobian(1,1) = 2*x1;
    Jacobian(1,2) = 1;
    Jacobian(2,1) = 1;
    Jacobian(2,2) = -2*x2;
    Jacobian(3,1) = 1;
    Jacobian(3,2) = 1;
    Jacobian(3,3) = 1;
    Jac = Jacobian;
end

% Newton's method that solves the function
% INPUTS: J - the Jacobian, F - The function to solve, 
%         N - Max iteration, ini - initial guess
function X = Newton(J,F,N,ini)
    x_old = ini;
    for i = 1:N
        Jac = J(x_old(1), x_old(2), x_old(3));
        Y = F(x_old(1), x_old(2), x_old(3));
        x_new = x_old - Jac\Y;
        x_old = x_new;
    end
    X = x_new;
end