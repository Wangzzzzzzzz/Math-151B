% run a test on the functions below
disp(heun(0.1,0,1,1,@f))


% Use function the following function as a test
% the actual soltuion is y = exp(t)
function dydt = f(t,y)
    dydt = y^2*exp(-t);
end

% Heun's Method
% input h, a, b, alpha (initial condition), func
function y = heun(h,a,b,alpha,func)
    t = a;
    w = alpha;
    N = (b-a)/h;
    for i = 1:N
        K1 = h/3 * func(t,w);
        K2 = 2/3 * h * func(t + h/3, w + K1);
        K3 = 3 * func(t + 2/3*h, w + K2);
        w = w + h/4 * ( func(t,w) + K3 );
        t = a + i*h;
    end
    y = w;
end