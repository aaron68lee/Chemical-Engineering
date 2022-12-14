
clc; close all; clear all;

% define function and params

h = 0.1; % step size

% Problem 1
y0(1, :) = [1, 1]; % x, x', x'' 
x0 = 0;
endT = 10;

% Problem 5
y0 = 1; % initial y
x0 = 0;
target = [1, 0; 0, 1];
endT = 1;
nodes = 10;


%[t, y] = RK(h, x0, y0, endT);
%[t, y] = shooting(h, x0, y0, target, endT);
[t, y] = FD(x0, y0, endT, target, nodes);

disp("Y Iterations:");
disp(y);
fprintf("Final y: \n");
disp(y(end));

%plot(x, y);
%plot(t, y(:, 1), t, y(:, 2));
plot(t, y(1:nodes+2), t, y(nodes+3:2*(nodes+2)));
legend('Y1', 'Y2'); %, "y'", "y''")
grid on 

%% Functions

%% Boundary Differential Equations Q6 Function and Jacobians

% returns Jacobian matrix
function J = BVP_J(Y, h)

    n = size(Y,1)/2;
    J = zeros(2*(n-2),2*(n-2));
    for i = 2:n-1
        ii = i+n;
        j = i-1;
        jj = n-2+j;
        J(j,j)=-2/h^2-2*Y(i);

        if j<n-2
            J(j,j+1)=1/h^2;
        end

        if j>1
            J(j,j-1)=1/h^2;
        end

        J(j,jj) = 1;
        J(jj,jj) = -2/h^2-1;
        if j < n-2
            J(jj,jj+1) = 1/h^2;
        end
        if j > 1
            J(jj, jj-1) = 1/h^2;
        end
        J(jj,j) = 2*Y(i);
    end
end

function F = BVP_F(Y,h)
    n = size(Y,1)/2;
    
    j = 0;
    F = zeros(2*(n-2),1);
    for i = 2:n-1
        ii = i+n;
        j = i-1;
        jj = n-2+j;
        F(j) = (Y(i+1)-2*Y(i)+Y(i-1)) / h^2-Y(i)^2+Y(ii);
        F(jj) = (Y(ii+1)-2*Y(ii)+Y(ii-1))/h^2+Y(i)^2-Y(ii);
    end
end

%% Finite Difference

function [x, y] = FD(x0, y0, endX, endY, nodes)

    %% Create matrix of equations
    %{
    A = zeros(nodes + 2, nodes + 2);
    A(1, 1) = 1;
    A(nodes + 2, nodes + 2) = 1;

    for i = 2:nodes
        A(i, [i-1, i, i+1]) = [1, -2, 1];
    end

    % Insert Boundary Conditions

    B = zeros(nodes + 2, 1);
    B(1) = 1;
    B(nodes + 2) = 0.5;
    B(2:nodes+1) = h^2 .* (2 .* x.^3);

    sol = A\B;
    %}

    %% Central Finite Difference Question 5
    %{
    x = linspace(x0, endX, nodes); % number line
    A = zeros(nodes, nodes);
    B = zeros(nodes, 1);

    % Boundary Conditions

    A(1, 1) = 1;
    A(nodes, nodes) = 1;
    B(1) = x0;
    B(nodes) = endY;

    for i = 2:nodes-1
        A(i, i-1) = h^2;
        A(i, i) = h^2;
        A(i, i+1) = h^2;
        B(i) = h^2;
    end
    %}

    %{
    J = zeros(nodes, nodes);

    nodes = nodes + 2;
    h = (endX - x0) / (nodes - 1);

    y = zeros(1, nodes);
    dy = zeros(1, nodes);
    d2y = zeros(1, nodes);
    
    x = linspace(x0, endX, nodes); % number line

    ERROR = 1e-5;
    err = 1;
    
    j = 1;
    
 
    while err > ERROR
        J(1,1) = 1;
        J(nodes, nodes) = 1;
        F(1) = y(j,1) - y0;
        F(nodes) = y(j, nodes) - endY;
        for i = 2:nodes-1
            J(i,i-1) = 1/h^2;
            J(i,i) = -2/h^2 - 6*y(j,i)^2;
            J(i,i+1) = 1/h^2;
            F(i) = (y(j,i+1) - 2*y(j,i) + y(j,i-1))/h^2 - 2*y(j,i)^3;
        end
        y(j+1,:) = y(j,:)-(inv(J)*F')';
        err=0;
        for i = 2:nodes-1
            err = err + (y(j+1,i)-y(j,i))^2;
        end
        j=j+1;
    end

    %}

    %% Central Finite Difference Question 6
    
   N = nodes + 2;
   Y = zeros(2*N,1);
   Y(1) = endY(1, 1);
   Y(N) = endY(1, 2);
   Y(N+1) = endY(2, 1);
   Y(2*N) = endY(2, 2);
   h = 1/(N-1);
   X = 0:h:1;
   err = 1;
   ERROR = 1e-5;
   inner=[(2:N-1)';(N+2 : 2*N-1)'];

   while(err > ERROR)
       Y_Old = Y;
       Y(inner) = Y(inner) - inv(BVP_J(Y,h)) * BVP_F(Y,h);
       err = sum((Y_Old-Y).^2);
   end

   x = X;
   y = Y;
    %% syms f1 f2

    f1 = @(y1, y2, y3) (y3 - 2*y2 + y1) / (h^2); % node at y2
    f2 = @(y1, y2, y3) (y3 - y1) / (2*h); % node at y2

end

%% Evaluate multiple function values for Systems
function [ret] = f(t, z1, z2)

    f1 = @(t, z1, z2) (z2);
    f2 = @(t, z1, z2) (-z1*z1 + t + z1 - z2); % equal to z2'

    ret = [f1(t, z1, z2), f2(t, z1, z2)];
end

%% Diff eqs for question 5
function dy = shoot(y)
    dy(1) = y(2);
    dy(2) = 2 * y(1)^3;
end

%% Shooting Method

function [x, y] = shooting(h, x0, y0, target, endX)

    %{
    ERROR = 1e-5;
    guess1 = 5; % for z1'(0)
    guess2 = 10;
    error = 10;
    x(1) = x0(1, 1);
    y = y0;

    [xx, zz] = RK(h, x0, [y0, guess1], endX);
    sol1 = zz(end, 1);

    [xx, zz] = RK(h, x0, [y0, guess2], endX);
    sol2 = zz(end, 1);

    guess(1) = guess1;
    guess(2) = guess2;
    sol(1) = sol1;
    sol(2) = sol2;

    i = 3;

    while error > ERROR
        
        % now use linear interpolation to find a new guess
        m = (guess(i-2) - guess(i-1)) / (sol(i-2) - sol(i-1));
        guess(i) = guess(i-1) + m * (target - sol(i-1));

        % Propogate Again with RK4
        [xx, zz] = RK(h, x0, [y0, guess(i)], endX);
        sol(i) = zz(end, 1);

        % Update Guesses and Error
        error = abs(sol(i) - target);
        i = i + 1;
    end
    x = xx;
    y = zz;
    %}
        

    nodes = 11;
    h = x0 / (nodes - 1);
    
    x = 1:h:2;
    dy = [0,-2]; % initial y' guesses
    
    j = 1;
    ERROR = 1e-5;
    err = 1;

    while err > ERROR % Euler propogation
        if (j > 2)
            dy(j) = dy(j-2) + (target - y_hit(j-2)) / (y_hit(j-1) - y_hit(j-2))*(dy(j-1) - dy(j-2));
        end
        y(1,:) = [1, dy(j)];
        for i=1:nodes-1
            y(i+1,:) = y(i,:)+ h * shoot(y(i,:));
        end
        y_hit(j) = y(end,1);
        err = abs(y_hit(j)-target);
        j=j+1; 
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
    Z(1, :) = y0; % initial concentrations

    % system of differential equations
    % Problem 1
    f1 = @(t, z1, z2) (z2);
    f2 = @(t, z1, z2) (-z1*z1 + t + z1 - z2); % equal to z2'

    % Problem 2
    f1 = @(t, z1, z2) (z2);
    f2 = @(t, z1, z2) (-2*(z1^3)); % equal to z2'

    for i = 1:steps

        % 4 groupings since RK4

        k1_1 = f1(t(i), Z(i, 1), Z(i, 2)); 
        k1_2 = f2(t(i), Z(i, 1), Z(i, 2));

        k2_1 = f1(t(i) + h/2, Z(i, 1) + k1_1*h/2, Z(i, 2) + k1_2*h/2);
        k2_2 = f2(t(i) + h/2, Z(i, 1) + k1_1*h/2, Z(i, 2) + k1_2*h/2);

        k3_1 = f1(t(i) + h/2, Z(i, 1) + k2_1*h/2, Z(i, 2) + k2_2*h/2);
        k3_2 = f2(t(i) + h/2, Z(i, 1) + k2_1*h/2, Z(i, 2) + k2_2*h/2);

        k4_1 = f1(t(i) + h, Z(i, 1) + k3_1*h, Z(i, 2) + k3_2*h);
        k4_2 = f2(t(i) + h, Z(i, 1) + k3_1*h, Z(i, 2) + k3_2*h);

        % Update X
        x(i+1) = x(i) + h;
        t(i+1) = t(i) + h;

        % Update Y
        Z(i+1, 1) = Z(i, 1) + h/6 * (k1_1 + 2*k2_1 + 2*k3_1 + k4_1);
        Z(i+1, 2) = Z(i, 2) + h/6 * (k1_2 + 2*k2_2 + 2*k3_2 + k4_2);

    end

    y = Z; % return final concentrations
end