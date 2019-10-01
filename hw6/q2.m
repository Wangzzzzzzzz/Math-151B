% MATH 151b, HW 6
% Question 2
x = [1;1;1;1;1;1;1];
disp(broyden(x,@f,7,10^-5))

function y=f(w)
N = size(w,1);
h = 1/(N+1);
result = zeros(N,1);
result(1) = 2*h^2*w(1)^3 + 2*w(1) - w(2) - 1/2;
result(N) = 2*h^2*w(N)^3 + 2*w(N) - w(N-1) - 1/3;
for i=2:(N-1)
    result(i) = 2*h^2*w(i)^3 + 2*w(i) - w(i-1) - w(i+1);
end
y = result;
end

% Broyden?s Method
function [xv,it]=broyden(x,f,n,tol)
% Broyden's method for solving a system of n non-linear equations
% in n variables.
%
% Example call: [xv,it]=broyden(x,f,n,tol)
% Requires an initial approximation column vector x. tol is required
% accuracy. User must define function f
% xv is the solution vector, parameter it is number of iterations
% taken. WARNING. Method may fail, for example, if initial estimates
% are poor.
%
fr=zeros(n,1); it=0; xv=x;
%Set initial Br
h = 1/(n+1);
Br=zeros(n);
Br(1,1) = 6*h^2*x(1)^2+2;
Br(1,2) = -1;
Br(n,n-1) = -1;
Br(n,n) = 6*h^2*x(n)^2+2;
for i=2:(n-1)
    Br(i,i-1)=-1;
    Br(i,i)= 6*h^2*x(i)^2+2;
    Br(i,i+1) = -1;
end
fr=feval(f, xv);
while norm(fr)>tol
  it=it+1;
  pr=-Br*fr;
  tau=1;
  xv1=xv+tau*pr; xv=xv1;
  oldfr=fr; fr=feval(f,xv);
  %Update approximation to Jacobian using Broydens formula
  y=fr-oldfr; oldBr=Br;
  oyp=oldBr*y-pr; pB=pr'*oldBr;
  for i=1:n
    for j=1:n
      M(i,j)=oyp(i)*pB(j);
    end
  end
  Br=oldBr-M./(pr'*oldBr*y);
end
end