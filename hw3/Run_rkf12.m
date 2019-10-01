% Math 151b
% Homework 3, Question 2

% run the program with TOL = 10^-4
disp(run_rkf12(0, 1, 1, 10^-4, 0.5, 10^-7, @f))

% a demo function, a = 0, b = 1, alpha  = 1
function dydt = f(t,y)
    dydt = y^2*exp(-t);
end

% Implement RKF12 method with Euler's method and Modified Euler Method
% INPUTS: 
% a,b - endpoints; alpha - initial condition; TOL - tolerance
% hmax - maximum step size; hmin - min step size
% func - function to be solved
function y = run_rkf12(a, b, alpha, TOL, hmax, hmin, func)
    t = a;
    w = alpha;
    h = hmax;
    FLAG = 1;
    disp([t, w])
    while FLAG == 1
        K1 = h * func(t,w);
        K2 = func(t + h, w + K1);
        K3 = h/2 * ( func(t,w) + K2 );
        % Euler: w_i+1 = w_i + K1
        % M_Euler: w_i+1 = w_i + K3
        R = 1/h * abs(K3 - K1);
        if R <= TOL
            t = t + h;
            w = w + K1;
            disp([t,w,h]);
        end
        q = (1/2) * (TOL/R);
        % Adjust the step size
        if q <= 0.1
            h = 0.1 * h; % prevent delta to be too small and h goes to 0
        elseif q >= 4
            h = 4 * h;  % prevent delta become to large and h increase to fast
        else
            h = q * h; % normal case, just set h = qh
        end
        % bound h by hmax
        if h > hmax
            h = hmax; 
        end
        % terminating conditions
        if t >= b
            FLAG = 0; % reach the end point
        elseif t + h > b
            h = b - t; % adjust the final step
        elseif h < hmin
            FLAG = 0; 
            fprintf("Minimun h exceeded, Procedure completed unsuccessfully.")
        end
    end
    y = w;
end
        
    
    