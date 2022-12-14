%% Name: Aaron Lee
%{
UID: 505 540 473
Date: 10/7/22
%}

clc; close all; clear all;

% test example
%{
A = [[1, 2, 1, -1];
     [3, 2, 4, 4];
     [4, 4, 3, 4];
     [2, 0, 1, 5]];
B = [5; 16; 22; 15];
%}

% homework example
A = [[-1, 1, 0, 0];
     [0, 1, -4, 1];
     [1, -4, 1, 0];
     [0, 0, 1, -1]];
B = [1; -6; -4; -1];

n = length(B);

if (det(A) ~= 0)
    %disp("Correct: ");
    %disp(linsolve(A, B));
    A = [A B]; % augmented matrix
    x = gauss(A, B); % Use Gaussian Elimination Method
    %x = LU(A(:, 1:4), B); % Use LU Decomposition Method
    disp("Program Solution: ");
    disp(x);
else
    disp("No solution");
end

%% Functions

% Gaussian Elimination
function x = gauss(A, B)
    n = length(B);
    x = zeros(n, 1);
    m = zeros(n, n);
    % Forward Elimination
    for i = 1:n-1 % cols
        for j = i+1:n % rows
            m(j, i) = A(j, i) / A(i, i);
            for k = i:n % i+1
                A(j, k) = A(j, k) - m(j, i)*A(i, k);
            end
            B(j) = B(j) - m(j, i)*B(i);
            %{
            disp(A);
            disp(B);
            fprintf("=============\n");
            %}
            %m = A(j, i) / A(i, i); % scale factor
            %A(j, :) = A(j, :) - m*A(i, :);
        end
        A(:, n+1) = B;
    end
    %x(n) = A(n, n+1) / A(n, n);

    
    % Back Substitution
    U = A(:, 1:n); % given upper triangular
    x = backSub(U, B);
    
end

% LU Decomposition
function x = LU(A, B)
    n = length(B);
    [L, U] = lu(A);

    for i = 1:n

       % Finding L
       for k = 1:i-1
           L(i,k)=A(i,k);
           for j = 1:k-1
               L(i,k) = L(i,k) - L(i,j)*U(j,k);
           end
           L(i,k) = L(i,k) / U(k,k);
       end
       
       % Finding U
       for k = i:n
           U(i,k) = A(i,k);
           for j = 1:i-1
               U(i,k) = U(i,k)-L(i,j)*U(j,k);
           end
       end
    end

    % Set diagonal of L to 1
    for i = 1:n
        L(i,i)=1;
    end

    %disp(L);
    %disp(U);
    y = forwardSub(L, B);
    x = backSub(U, y);
end

% Backwards Substitution
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

% Forward Substitution
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