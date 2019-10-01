% Math 151b
% Homework 3, Question 2
time_rkf12 = zeros(1,7);
time_rkf45 = zeros(1,7);

for i=1:7
    tic
    disp(run_rkf12(0, 1, 1, 10^-i, 0.5, 10^-10, @f))
    time_rkf12(i) = toc;

    tic
    disp(run_rkf45(0, 1, 1, 10^-i, 0.5, 10^-10, @f))
    time_rkf45(i) = toc;
end
% make the plot
figure;
plot((1:7),time_rkf12,'Linewidth', 1.1);
hold on;
plot((1:7),time_rkf45,'Linewidth', 1.1);
ylim([0,0.1]);
xlabel('accuracy in 10^{-x}');
ylabel('Running Time (s)');
legend({'RKF12','RKF45'},'Location','northwest')
title('Plot of Performance against Accuracy');
grid on;
hold off;

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
    %disp([t, w])
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
            %disp([t,w,h]);
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
     

% Implement RKF45 method
% INPUTS: 
% a,b - endpoints; alpha - initial condition; TOL - tolerance
% hmax - maximum step size; hmin - min step size
% func - function to be solved
function y = run_rkf45(a, b, alpha, TOL, hmax, hmin, func)
    t = a;
    w = alpha;
    h = hmax;
    FLAG = 1;
    %disp([t, w])
    while FLAG == 1
        K1 = h * func(t,w);
        K2 = h * func(t + h/4, w + K1/4);
        K3 = h * func(t + 3/8*h, w + 3/32*K1 + 9/32*K2 );
        K4 = h * func(t + 12/13*h, w + 1932/2197*K1 - 7200/2197*K2 + 7296/2197*K3);
        K5 = h * func(t + h, w + 439/216*K1 -8*K2 + 3680/513*K3 - 845/4104*K4);
        K6 = h * func(t + h/2, w - 8/27*K1 + 2*K2 - 3544/2565*K3 + 1859/4104*K4 - 11/40*K5);
        R = 1/h * abs(1/360*K1 - 128/4275*K3 - 2197/75240*K4 + 1/50*K5 + 2/55*K6);
        if R <= TOL
            t = t + h;
            w = w + 25/216*K1 + 1408/2565*K3 + 2197/4104*K4 - 1/5*K5;
            %disp([t,w,h]);
        end
        q = 0.84 * (TOL/R)^(1/4);
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
    