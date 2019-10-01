% collect data
hs = (2:50).^-1;
i = 1;
err_heun = zeros(1,size(hs,1));
err_modi_euler = zeros(1,size(hs,1));
err_euler = zeros(1,size(hs,1));
for h = hs
    err_heun(i) = abs(heun(h,0,1,1,@f)-sol(1));
    err_modi_euler(i) = abs(modi_euler(h,0,1,1,@f)-sol(1));
    err_euler(i) = abs(euler(h,0,1,1,@f)-sol(1));
    i = i+1;
end
% make the plot
figure;
plot(hs,err_heun,'Linewidth', 1.1);
hold on;
plot(hs,err_modi_euler,'Linewidth', 1.1);
plot(hs,err_euler,'Linewidth', 1.1);
xlabel('h');
ylabel('Error');
legend({'Heun''s','Modified Euler', 'Euler'},'Location','northwest')
title('Plot of Error against h');
grid on;
hold off;

% Use function the following function as a test
% the actual soltuion is y = exp(t)
function dydt = f(t,y)
    dydt = y^2*exp(-t);
end

% solution
function s = sol(t)
    s = exp(t);
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

% Modified Euler's Method
% input h, a, b, alpha (initial condition), func
function y = modi_euler(h,a,b,alpha,func)
    t = a;
    w = alpha;
    N = (b-a)/h;
    for i = 1:N
        K1 = h * func(t,w);
        K2 = func(t + h, w + K1);
        w = w + h/2 * ( func(t,w) + K2 );
        t = a + i*h;
    end
    y = w;
end

% Euler's Method
% input h, a, b, alpha (initial condition), func
function y = euler(h,a,b,alpha,func)
    t = a;
    w = alpha;
    N = (b-a)/h;
    for i = 1:N
        w = w + h * func(t,w);
        t = a + i*h;
    end
    y = w;
end