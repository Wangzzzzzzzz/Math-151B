% Math 151b, Homework 4, Question 2(b),2(c)
% Wang, Zheng 

% test with the function given in part (a)
% the correct solution is about 3.16177
pred_cor(0.1, 0, 1, 0, 0, @dudt, @dvdt)

% du/dt
function f_u = dudt(t,u,v)
    f_u = t*0 + u*0 + v;
end

% dv/dt
function f_v = dvdt(t,u,v)
    f_v = 4*u + 6*exp(-t) + 0*v;
end

% function of predictor-corrector method
% input h, a, b, alpha_u (initial condition of u), alpha_v, f_u, f_v
function y = pred_cor(h,a,b,alpha_u,alpha_v,f_u,f_v)
    t = a;
    U = alpha_u;
    V = alpha_v;
    N = (b-a)/h;
    for i = 1:N
        % Predictor Step
        Ku_1 = U + h/2 * f_u(t,U,V);
        Kv_1 = V + h/2 * f_v(t,U,V);
        Ku_2 = f_u(t + h/2, Ku_1, Kv_1);
        Kv_2 = f_v(t + h/2, Ku_1, Kv_1);
        U_temp = U; % store U_i
        V_temp = V; % store V_i
        t_temp = t; % store t_i
        U = U + h*Ku_2; % update to U_i+1 (prediction)
        V = V + h*Kv_2; % update to V_i+1 (prediction)
        t = a + i*h;    % update to t_i+1
        % Corrector step
        % Correct U and V with the one-step implicit method and pass 
        % to next iteration
        U = U_temp + h/2 * f_u(t,U,V) + h/2 * f_u(t_temp,U_temp,V_temp);
        V = V_temp + h/2 * f_v(t,U,V) + h/2 * f_v(t_temp,U_temp,V_temp);
    end
    y = U;
end