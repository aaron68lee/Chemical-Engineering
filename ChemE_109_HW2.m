%% Name: Aaron Lee
%{
UID: 505 540 473
Date: 10/26/22
%}

clc; close all; clear all;


% homework question 2: LU

A = [[0, 0, 1, -1];
     [0, 1, -4, 1];
     [1, -4, 1, 0];
     [-1, 1, 0, 0]];
B = [7; -18; 12; 3];


% homework question 4: Gauss-Seidel and Jacobi
%{
A = [[5, 0.1, 0.3, -1];
     [0.5, 3, -0.4, 0.1];
     [0.1, -0.4, 5, 0];
     [-0.2, 0.1, 0, 5]];
B = [1; 2; 3; 4];
%}

%{
A = [[5, -0.3, 0.1];
     [0.2, -2, -1];
     [0.3, -0.1, -4]];
B = [1; 2; 3];
%}

x0 = [0; 0; 0; 0];
n = length(B);

[L, U, P] = lu(A);

if (det(A) ~= 0)
    %disp("Correct: ");
    %disp(linsolve(A, B));
    
    %[P, x] = LUSolver(A, B);
    %x = P*x;
    %x = JacobiSolver(A, B, x0); % Use Jacobian Method
    x = GaussSeidelSolver(A, B, x0); % Use Gauss-Seidel Method
    
    disp("Program Solution: ");
    disp(x);
else
    disp("No solution");
end

%% =================== Functions ===================

%% Gaussian Elimination
function x = gauss(A, B)
    n = length(B);
    x = zeros(n, 1);
    m = zeros(n, n);

    % Forward Elimination
    for i = 1:n-1 % cols
        for j = i+1:n % rows
            m(j, i) = A(j, i) / A(i, i); % find row ratio
            for k = i:n % i+1
                A(j, k) = A(j, k) - m(j, i)*A(i, k);
            end
            % modify B vector
            B(j) = B(j) - m(j, i)*B(i);
        end
        A(:, n+1) = B; % replace last row of A with B
    end
    %x(n) = A(n, n+1) / A(n, n);

    
    % Back Substitution
    U = A(:, 1:n); % given upper triangular
    x = backSub(U, B);
    
end

%% LU Decomposition
function [P, x] = LUSolver(A, B)
    n = length(B);
    
    L = zeros(n, n);
    U = zeros(n, n);
    P = eye(n, n);
    [L, U, P] = lu(A);

    % Partial Pivoting
    A = A*P;
    
    %{
    for i = 1:n-1
        [M, I] = max(abs(A(i, :)));
        I = I + i - 2;
        if I ~= i % Exchange rows i and k, ignoring columns 1 through i-1 in each row.
            temp = A(i, :);
            A(i, :) = A(I, :);
            A(I, :) = temp;
        end
    end
    %}

    for i = 1:n
       % Finding L
       for k = 1:i-1
           L(i,k) = A(i,k);
           for j = 1:k-1
               L(i,k) = L(i,k) - L(i,j)*U(j,k);
           end
           L(i,k) = L(i,k) / U(k,k);
       end
       
       % Finding U
       for k = i:n
           U(i,k) = A(i,k);
           for j = 1:i-1
               U(i,k) = U(i,k) - L(i,j)*U(j,k);
           end
       end
    end

    % Set diagonal of L to 1
    for i = 1:n
        L(i,i)=1;
    end

    disp(L);
    disp(U);
    y = forwardSub(L, B);
    x = backSub(U, y);
end


%% Gauss-Seidel Method
function x = GaussSeidelSolver(A, b, x0)
    error = 0.5;
    plotGauss = [];
    count = 0;
    n = length(b);
    normVal = 10;
    x = x0;

    %x = rand(5, 1);
    while normVal > error
        x_old = x;
        
        % print results of current iteration
        fprintf("Iteration " + num2str(count) + ":\n");
        disp(x);

        for i=1:n

            sigma=0;
            
            for j=1:i-1
                    sigma = sigma + A(i,j)*x(j);
            end

            for j=i+1:n
                    sigma = sigma + A(i,j)*x_old(j);
            end
            
            x(i) = (1/A(i,i)) * (b(i)-sigma);
        end
        
        count = count + 1;
        normVal = norm(x_old - x);
        plotGauss = [plotGauss;normVal];
    end

    %{
    hold on
    plot(1:5:count, plotGauss(1:5:count),'LineWidth',2)
    text(count, 0.2, '\downarrow')
    text(count,0.3,'Gauss Seidel')
    legend('Gauss Seidel Method')
    ylabel('Error Value')
    xlabel('Number of iterations')
    title('Iterative Methods')
    hold off
    %}

end

%% Jacobi Method
function x = JacobiSolver(A, b, x0)

    error = 0.5;
    plotJacobi = [];
    count = 0;
    n = length(b);
    normVal = 10;
    x = x0;

    while normVal > error
        % log previous iteration in temp storage
        x_old = x;

        % print current iteration
        fprintf("Iteration " + num2str(count) + ":\n");
        disp(x_old);

        for i = 1:n
            sigma=0;
            
            for j = 1:n
                if j ~= i
                    sigma = sigma + A(i,j)*x(j);
                end
            end
            
            x(i) = (1/A(i,i)) * (b(i)-sigma);
        end
        
        count = count + 1;
        normVal = norm(x_old-x);
        plotJacobi = [plotJacobi;normVal];
    end

    %{
    hold on
    plot(1:5:count, count(1:5:count),'LineWidth',2)
    text(count,0.3,'\downarrow')
    text(count,0.4,'Jacobi Method')
    legend('Jacobi Method')
    ylabel('Error Value')
    xlabel('Number of iterations')
    title('Iterative Methods')
    hold off
    %}

end

%% Backwards Substitution
function x = backSub(U, B)
    n = length(B)
    x = zeros(n, 1);
    for i = n:-1:1
        temp = B(i);
        for j = (i+1):n
            temp = temp - x(j) * U(i, j);
        end
        x(i) = temp / U(i, i);
    end
end

%% Forward Substitution
function x = forwardSub(L, B)
    n = length(B);
    x = zeros(n, 1);
    for i = 1:n
        temp = B(i);
        for j = 1:i-1
            temp = temp - x(j) * L(i, j);
        end
        x(i) = temp / L(i, i);
    end
end