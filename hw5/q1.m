% Math 151b HW5
% Question 1

% generate sequence of the real and imaginary axis
re = -3:0.01:3;
img = -3:0.01:3;

% run contour map, plot only Z from 0 to 1, shown in blue
figure;
contourf(re,img,Q_euler(re,img),[0 1]);
title("Absolute Stability Region of Euler's Method");
xlabel("Real Line")
ylabel("Imaginary Line")
grid on;
pbaspect([1 1 1]);


% run contour map, plot only Z from 0 to 1, shown in blue
figure;
contourf(re,img,Q_mid(re,img),[0 1]);
title("Absolute Stability Region of Midpoint Method");
xlabel("Real Line")
ylabel("Imaginary Line")
grid on;
pbaspect([1 1 1]);


% Generate a grid of |Q(lambda h)| 
% otherwise
% real: real sequence, img: imaginary sequence
function z = Q_euler(rel,img)
    M = zeros(size(img,2), size(rel,2));
    for j = 1:size(rel,2)
        for k = 1:size(img,2)
            lambda_h = rel(j)+img(k)*1i;
            M(k,j) = abs(1+lambda_h);
        end
    end
    z = M;
end
    
    
% Generate a grid of |Q(lambda h)| 
% otherwise
% real: real sequence, img: imaginary sequence
function z = Q_mid(rel,img)
    M = zeros(size(img,2), size(rel,2));
    for j = 1:size(rel,2)
        for k = 1:size(img,2)
            lambda_h = rel(j)+img(k)*1i;
            M(k,j) = abs(1+lambda_h+0.5*lambda_h^2);
        end
    end
    z = M;
end
    
    