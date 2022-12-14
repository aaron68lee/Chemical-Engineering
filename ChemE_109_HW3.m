%% Name: Aaron Lee
%{
UID: 505 540 473
Date: 11/04/22

Uncomment certain sections of code to execute different parts of the HW
%}

clc; close all; clear all;

% define function and params

h = 0.1; % step size

x0(1, :) = [0];
y0(1, :) = [1 0 0]; % initial guesses for sys differential equations

endX = 10;


%[x, y] = RK(h, x0, y0, endX);
%[x, y] = EPC(h, x0, y0, endX);



x0 = 2;
y0 = 2;
endX = 10;
h = 0.01;
[x, y] = EE2(h, x0, y0, endX);

disp("Y Iterations:");
disp(y);
fprintf("Final y: \n");
disp(y(end, :));

plot(x, y);
%plot(x, y(:, 1), x, y(:, 2));
%plot(x, y(:, 1), x, y(:, 2), x, y(:, 3));
legend('First:', "Second:", "Third:")
grid on

%% =================== Functions ===================

%% Evaluate multiple function values for Systems
function [ret] = f(t, x, y, z)
    
    % part 1
    f0 = @(x, y) (x + y)/2; 
    
    % part 2
    %f1 = @(y1, y2) (-0.5*y1); 
    %f2 = @(y1, y2) (4 - 0.3*y2 - 0.1*y1);

    % part 3

    k1 = 1;
    k2 = 2;
    k3 = 3;

    f1 = @(t, a, b, c) (-k1*a);
    f2 = @(t, a, b, c) (k1*a - k2*b + k3*c);
    f3 = @(t, a, b, c) (k2*b - k3*c);
    
    %ret = [f1(t, x, y, z), f2(t, x, y, z), f3(t, x, y, z)]; %f0(x, y)
end

%% EE Problem 6
function [x, y] = EE2(h, x0, y0, endX)
    
    steps = ceil((endX - x0) / h);
    x(1) = x0;
    y(1) = y0;

    f0 = @(x, y) (x + y)/2; 

    for i = 1:steps
        x(i+1) = x(i) + h;
        fun = f0(x(i), y(i));
        y(i+1) = y(i) + h * fun;
    end
end

%% Explicit Euler's Method
function [x, y] = EE(h, x0, y0, endX)
    
    steps = ceil((endX - x0(1, 1)) / h);
    x(1) = x0(1, 1);
    y(1, :) = y0;

    %f0 = @(x, y) (x + y)/2;

    for i = 1:steps
        x(i+1) = x(i) + h;
        fun = f0(x(i), y(i, :));
        y(i+1, :) = y(i, :) + h * fun(2:end);
    end
end

%% Euler Predictor-corrector Method
function [t, y] = EPC(h, x0, y0, endX)
    
    steps = ceil((endX - x0(1, 1)) / h);

    t(1) = x0(1, 1);
    y(1, :) = y0;

    k1 = 1;
    k2 = 2;
    k3 = 3;

    f1 = @(a, b, c) (-k1*a);
    f2 = @(a, b, c) (k1*a - k2*b + k3*c);
    f3 = @(a, b, c) (k2*b - k3*c);

    for i = 1:steps
        t(i+1) = t(i) + h;
        yP(i+1, :) = y(i, :) + h * f(t(i), y(i, 1), y(i, 2), y(i, 3));
        fun = f(t(i), y(i, 1), y(i, 2), y(i, 3));
        fun2 = f(t(i+1), yP(i+1, 1), yP(i+1, 2), yP(i+1, 3));
        y(i+1, :) = y(i, :) + h/2 * (fun + fun2);
    end
end

%% Newton's Method

function x = newton(h, x0, y0, endX)

    counter = 1;
    MAX_ITER = 1000;

    steps = ceil((endX - x0(1, 1)) / h);
    x(1) = x0(1, 1);
    y(1, :) = y0;

    tempV = x(end);

    f1 = @(y1, y2) (-0.5*y1);
    df1 = @(y1, y2) (-0.5*y1); % derivative of function

    % Newton's Method
    while abs(F) > ERROR
        counter = counter + 1;
        F = f1(x, y);
        dF = df1(x, y);
        x(end+1) = x(end) - F / dF;
        if (counter > MAX_ITER)
            fprintf("\nMax iterations reached without convergence.\n");
            break;
        end
    end
end

%% 4th Order Runge-Kutta
function [x, y] = RK(h, x0, y0, endX)
% params: step size, x initial, y initial, final X, function
% do fourth order approximation here
    
    steps = ceil((endX - x0) / h);
    x(1) = x0;
    y(1, :) = y0;

    t(1) = x0;
    C(1, :) = y0; % initial concentrations

    % system of differential equations
    
    k1 = 1;
    k2 = 2;
    k3 = 3;

    f1 = @(t, a, b, c) (-k1*a);
    f2 = @(t, a, b, c) (k1*a - k2*b + k3*c);
    f3 = @(t, a, b, c) (k2*b - k3*c);

    for i = 1:steps
        
        %{
        k1_1 = f1(x(i), y(i, 1)); % dt(trial) * D * Laplacian_2D(grid_temp, h, stencil);
        k2_1 = f1(x(i) + h/2, y(i, 1) + k1_1*h/2);
        k3_1 = f1(x(i) + h/2, y(i, 1) + k2_1*h/2);
        k4_1 = f1(x(i) + h, y(i, 1) + k3_1*h);
        
        k1_2 = f2(x(i), y(i, 2));
        k2_2 = f2(x(i) + h/2, y(i, 2) + k1_2*h/2);
        k3_2 = f2(x(i) + h/2, y(i, 2) + k2_2*h/2);
        k4_2 = f2(x(i) + h, y(i, 2) + k3_2*h);
        %}


        k1_1 = f1(t(i), C(i, 1), C(i, 2), C(i, 3)); % dt(trial) * D * Laplacian_2D(grid_temp, h, stencil);
        k1_2 = f2(t(i), C(i, 1), C(i, 2), C(i, 3));
        k1_3 = f3(t(i), C(i, 1), C(i, 2), C(i, 3));

        k2_1 = f1(t(i) + h/2, C(i, 1) + k1_1*h/2, C(i, 2) + k1_2*h/2, C(i, 3) + k1_3*h/2);
        k2_2 = f2(t(i) + h/2, C(i, 1) + k1_1*h/2, C(i, 2) + k1_2*h/2, C(i, 3) + k1_3*h/2);
        k2_3 = f3(t(i) + h/2, C(i, 1) + k1_1*h/2, C(i, 2) + k1_2*h/2, C(i, 3) + k1_3*h/2);

        k3_1 = f1(t(i) + h/2, C(i, 1) + k2_1*h/2, C(i, 2) + k2_2*h/2, C(i, 3) + k2_3*h/2);
        k3_2 = f2(t(i) + h/2, C(i, 1) + k2_1*h/2, C(i, 2) + k2_2*h/2, C(i, 3) + k2_3*h/2);
        k3_3 = f3(t(i) + h/2, C(i, 1) + k2_1*h/2, C(i, 2) + k2_2*h/2, C(i, 3) + k2_3*h/2);

        k4_1 = f1(t(i) + h, C(i, 1) + k3_1*h, C(i, 2) + k3_2*h, C(i, 3) + k3_3*h);
        k4_2 = f2(t(i) + h, C(i, 1) + k3_1*h, C(i, 2) + k3_2*h, C(i, 3) + k3_3*h);
        k4_3 = f3(t(i) + h, C(i, 1) + k3_1*h, C(i, 2) + k3_2*h, C(i, 3) + k3_3*h);

        % Update X
        x(i+1) = x(i) + h;
        t(i+1) = t(i) + h;

        % Update Y
        %y(i+1, 1) = y(i, 1) + h/6 * (k1_1 + 2*k2_1 + 2*k3_1 + k4_1);
        %y(i+1, 2) = y(i, 2) + h/6 * (k1_2 + 2*k2_2 + 2*k3_2 + k4_2);
        C(i+1, 1) = C(i, 1) + h/6 * (k1_1 + 2*k2_1 + 2*k3_1 + k4_1);
        C(i+1, 2) = C(i, 2) + h/6 * (k1_2 + 2*k2_2 + 2*k3_2 + k4_2);
        C(i+1, 3) = C(i, 3) + h/6 * (k1_3 + 2*k2_3 + 2*k3_3 + k4_3);

    end

    y = C; % return final concentrations
end
