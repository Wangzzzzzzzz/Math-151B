% Math 151b HW5
% Question 2

N = 100;
y = solve_bvp(N);
x = linspace(0,1,N);

% show the plot
plot(x,f(x));
hold on;
plot(x,y);
grid on;
title("Plot of solution to BVP");
xlabel("x")
ylabel("y")
legend('Estimate','Actual','Location','northwest')
hold off;

% display y(0)
disp(y(1))


% actual answer
function y = f(x)
    y = (2*exp(4)*x + 2*x + exp(4-2*x) - exp(2*x))/(2+2*exp(4));
end

% function that solve the BVP in the question 
% input: N - number of grids
function w = solve_bvp(N)
    h = 1/N; 
    % fill A
    A = zeros(N,N);
    A(1,1) = -3;
    A(1,2) = 4;
    A(1,3) = -1;
    for i = 2:(N-1)
        A(i,i-1) = 1;
        A(i,i) = -(2+4*h^2);
        A(i,i+1) = 1;
    end
    A(N,N-1) = 1;
    A(N,N) = -(2+4*h^2);
    % fill b
    b = zeros(N,1);
    b(1) = 0;
    for i = 2:(N-1)
        b(i) = -4*h^2*(i-1)*h;
    end
    b(N) = -4*h^2*(N-1)*h - 1;
    % solve for w
    w = A\b;
end
    
        
    
