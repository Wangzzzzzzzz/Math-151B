% MATH 151b, HW 6
% Question 3
x = [1;1;1];
fprintf('Use Tolerance = 0.01\n')
solution = Ste_Dec(@F,@J,@g,x,0.01,100000000);
fprintf('Solution find is:\n')
disp(solution)
fprintf('Check the solution is close to actual solution, F(x) is\n')
disp(F(solution))
fprintf('\n\n');
fprintf('Use Tolerance = 10^-5\n')
solution = Ste_Dec(@F,@J,@g,x,10^-5,100000000);
fprintf('Solution find is:\n')
disp(solution)
fprintf('Check the solution is close to actual solution, F(x) is\n')
disp(F(solution))


function Y = F(x)
    y = zeros(3,1);
    y(1) = x(1)^3 + x(1)^2*x(2) - x(1)*x(3) + 6;
    y(2) = exp(x(1)) + exp(x(2)) - x(3);
    y(3) = x(2)^2 - 2*x(1)*x(3) - 4;
    Y = y;
end

function Jacb = J(x)
    Jac = zeros(3);
    Jac(1,1) = 3*x(1)^2 + 2*x(1)*x(2) - x(3);
    Jac(1,2) = x(1)^2;
    Jac(1,3) = -x(1);
    Jac(2,1) = exp(x(1));
    Jac(2,2) = exp(x(2));
    Jac(2,3) = -1;
    Jac(3,1) = -2*x(3);
    Jac(3,2) = 2*x(2);
    Jac(3,3) = -2*x(1);
    Jacb = Jac;
end

function y = g(x)
    f1 = x(1)^3 + x(1)^2*x(2) - x(1)*x(3) + 6;
    f2 = exp(x(1)) + exp(x(2)) - x(3);
    f3 = x(2)^2 - 2*x(1)*x(3) - 4;
    y = f1^2 + f2^2 + f3^2;
end

function result = Ste_Dec(F,J,g,ini,tol,max_iter)
    x = ini;
    k = 1;
    while k <= max_iter
        g1 = g(x);
        z = 2*J(x).'*F(x);
        z0 = norm(z);
        if z0 == 0
            result = x;
            fprintf('Iteration number:');
            disp(k);
            return
        end
        z = z/z0;
        alpha1 = 0;
        alpha3 = 1;
        g3 = g(x-alpha3*z);
        while g3 >= g1
            alpha3 = alpha3/2;
            g3 = g(x-alpha3*z);
            if alpha3 < tol/2
                fprintf('No likely improvement\n');
                result = x;
                fprintf('Iteration number:');
                disp(k);
                return
            end   
        end
        alpha2 = alpha3/2;
        g2 = g(x-alpha2*z);
        % solve for minimum of the interpolation function
        h1 = (g2-g1)/alpha2;
        h2 = (g3-g2)/(alpha3-alpha2);
        h3 = (h2-h1)/alpha3;
        alpha0 = 0.5*(alpha2-h1/h3);
        g0 = g(x-alpha0*z);
        if g3 <= g0
            alpha = alpha3;
            g_val = g3;
        else
            alpha = alpha0;
            g_val = g0;
        end
        x = x-alpha*z;
        if abs(g_val-g1) < tol
            result = x;
            fprintf('Iteration number:');
            disp(k);
            return
        end
        k = k+1;
    end
    result = x;
    fprintf('Reach max iteration\n');
end
        
        
    