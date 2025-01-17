fprintf("h = 0.5, Euler's method approximation result is y(1)=%f, Taylor method of order 2 approximation result is y(1)=%f.\n", euler(0.5), taylor(0.5));
fprintf("h = 0.1, Euler's method approximation result is y(1)=%f, Taylor method of order 2 approximation result is y(1)=%f.\n", euler(0.1), taylor(0.1));
fprintf("h = 0.01, Euler's method approximation result is y(1)=%f, Taylor method of order 2 approximation result is y(1)=%f.\n", euler(0.01), taylor(0.01));

function y1 = euler(h)
    y = 1;
    ts = linspace(0,1,int32(1/h)+1);
    for t = ts(1:end-1)
        y = y + h * (y^2*exp(-t));
    end
    y1 = y;
end

function y1 = taylor(h)
    y = 1;
    ts = linspace(0,1,int32(1/h)+1);
    for t = ts(1:end-1)
        y = y + h * (y^2*exp(-t)) + h^2/2 * (-y^2*exp(-t) + 2*y^3*exp(-2*t));
    end
    y1=y;
end