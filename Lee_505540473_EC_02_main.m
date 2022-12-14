
%% Extra Credit Problem 2 - Newton's Method
% Aaron Lee
% UID: 505540473
% Objective: Use Newton's Method to find zeros of a complex function

clc; clear; close all;

deltas = [1e-6, 1e-8, 1e-10];
% create handle to polynomial function
f = @polynomial;

for delta = 1:3 % run through each value of delta
    for i = 1.43:0.01:1.71 % run through guess-root range [1.43, 1.71]
        [root, evals] = Newton(f, i, deltas(delta), 100);
        fprintf("x0 = %1.2f, evals = %i, xc = %2.6f\n", i, evals, root); %, deltas(delta));
    end
    fprintf("\n\n");
end


%% function implementations

function y = polynomial(x)
% a polynomial function to test out Newton's Method on
% accepts one parameter, x, and returns y

y = 816*x^3 - 3835*x^2 + 6000*x - 3125;

end

function [xc, fEvals] = Newton(f, x0, delta, fEvalMax)
% f is a HANDLE to a continuous function, f(x), of a single variable.
% x0 is an initial guess to a root of f.
% delta is a positive real number.
% fEvalsMax is a positive integer >= 2 that indicates maximum
% number of f-evaluation allowed.
%
% Newton's method is repeatedly applied until the current iterate, xc,
% has the property that |f(xc)| <= delta. If that is not the case
% after fEvalsMax function evaluations, then xc is the current iterate.
%
% This routine computes the derivative of f at each iterate xc by 1
% using a central difference approximation with small perturbation size.
%
% fEvals is the number of f-evaluations required to obtain xc.
    
    % variable initilization 

    count = 0;
    xc = x0;
    h = 1e-6;
    
    while count < fEvalMax && abs(f(xc)) > delta 
        % keep iterating Newton's method until close to actual root or max evals exceeded
        
        df = (f(xc + h) - f(xc - h)) / (2*h); % calculate f'(xc)
        x1 = xc - f(xc) / df;
        xc = x1; % update new root guess value
        count = count + 1;
    end
    
    fEvals = count; % return eval count
end

%% Discussion Questions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
comment start

1) Console Output where max evaluations is set to 100
    
    Delta = 1e-6:

        x0 = 1.43, evals = 5, xc = 1.470588
        x0 = 1.44, evals = 4, xc = 1.470588
        x0 = 1.45, evals = 4, xc = 1.470588
        x0 = 1.46, evals = 3, xc = 1.470588
        x0 = 1.47, evals = 2, xc = 1.470588
        x0 = 1.48, evals = 3, xc = 1.470588
        x0 = 1.49, evals = 4, xc = 1.470588
        x0 = 1.50, evals = 6, xc = 1.470588
        x0 = 1.51, evals = 18, xc = 1.666667
        x0 = 1.52, evals = 8, xc = 1.470588
        x0 = 1.53, evals = 4, xc = 1.562500
        x0 = 1.54, evals = 3, xc = 1.562500
        x0 = 1.55, evals = 3, xc = 1.562500
        x0 = 1.56, evals = 2, xc = 1.562500
        x0 = 1.57, evals = 2, xc = 1.562500
        x0 = 1.58, evals = 3, xc = 1.562500
        x0 = 1.59, evals = 3, xc = 1.562500
        x0 = 1.60, evals = 4, xc = 1.562500
        x0 = 1.61, evals = 6, xc = 1.666667
        x0 = 1.62, evals = 8, xc = 1.470588
        x0 = 1.63, evals = 7, xc = 1.666667
        x0 = 1.64, evals = 5, xc = 1.666667
        x0 = 1.65, evals = 4, xc = 1.666667
        x0 = 1.66, evals = 3, xc = 1.666667
        x0 = 1.67, evals = 3, xc = 1.666667
        x0 = 1.68, evals = 3, xc = 1.666667
        x0 = 1.69, evals = 4, xc = 1.666667
        x0 = 1.70, evals = 4, xc = 1.666667
        x0 = 1.71, evals = 5, xc = 1.666667

    Delta = 1e-8:

        x0 = 1.43, evals = 5, xc = 1.470588
        x0 = 1.44, evals = 5, xc = 1.470588
        x0 = 1.45, evals = 4, xc = 1.470588
        x0 = 1.46, evals = 4, xc = 1.470588
        x0 = 1.47, evals = 2, xc = 1.470588
        x0 = 1.48, evals = 4, xc = 1.470588
        x0 = 1.49, evals = 5, xc = 1.470588
        x0 = 1.50, evals = 6, xc = 1.470588
        x0 = 1.51, evals = 18, xc = 1.666667
        x0 = 1.52, evals = 9, xc = 1.470588
        x0 = 1.53, evals = 4, xc = 1.562500
        x0 = 1.54, evals = 3, xc = 1.562500
        x0 = 1.55, evals = 3, xc = 1.562500
        x0 = 1.56, evals = 2, xc = 1.562500
        x0 = 1.57, evals = 2, xc = 1.562500
        x0 = 1.58, evals = 3, xc = 1.562500
        x0 = 1.59, evals = 4, xc = 1.562500
        x0 = 1.60, evals = 4, xc = 1.562500
        x0 = 1.61, evals = 6, xc = 1.666667
        x0 = 1.62, evals = 8, xc = 1.470588
        x0 = 1.63, evals = 7, xc = 1.666667
        x0 = 1.64, evals = 5, xc = 1.666667
        x0 = 1.65, evals = 4, xc = 1.666667
        x0 = 1.66, evals = 4, xc = 1.666667
        x0 = 1.67, evals = 3, xc = 1.666667
        x0 = 1.68, evals = 4, xc = 1.666667
        x0 = 1.69, evals = 4, xc = 1.666667
        x0 = 1.70, evals = 5, xc = 1.666667
        x0 = 1.71, evals = 5, xc = 1.666667


    Delta = 1e-10:

        x0 = 1.43, evals = 5, xc = 1.470588
        x0 = 1.44, evals = 5, xc = 1.470588
        x0 = 1.45, evals = 5, xc = 1.470588
        x0 = 1.46, evals = 4, xc = 1.470588
        x0 = 1.47, evals = 3, xc = 1.470588
        x0 = 1.48, evals = 4, xc = 1.470588
        x0 = 1.49, evals = 5, xc = 1.470588
        x0 = 1.50, evals = 7, xc = 1.470588
        x0 = 1.51, evals = 18, xc = 1.666667
        x0 = 1.52, evals = 9, xc = 1.470588
        x0 = 1.53, evals = 4, xc = 1.562500
        x0 = 1.54, evals = 4, xc = 1.562500
        x0 = 1.55, evals = 3, xc = 1.562500
        x0 = 1.56, evals = 3, xc = 1.562500
        x0 = 1.57, evals = 3, xc = 1.562500
        x0 = 1.58, evals = 3, xc = 1.562500
        x0 = 1.59, evals = 4, xc = 1.562500
        x0 = 1.60, evals = 4, xc = 1.562500
        x0 = 1.61, evals = 7, xc = 1.666667
        x0 = 1.62, evals = 9, xc = 1.470588
        x0 = 1.63, evals = 8, xc = 1.666667
        x0 = 1.64, evals = 6, xc = 1.666667
        x0 = 1.65, evals = 5, xc = 1.666667
        x0 = 1.66, evals = 4, xc = 1.666667
        x0 = 1.67, evals = 3, xc = 1.666667
        x0 = 1.68, evals = 4, xc = 1.666667
        x0 = 1.69, evals = 5, xc = 1.666667
        x0 = 1.70, evals = 5, xc = 1.666667
        x0 = 1.71, evals = 5, xc = 1.666667

2) Deltas

    Decreasing the size of delta generally increases the number of
    evaluations to reach the root. Looking at the last 5 entries in each
    set of delta values used, namely root guesses ranging from 1.67 to
    1.71, we see that the number of evaluations needed were as follows:
        
        Delta = 1e-6:        Evals needed: 3, 3, 4, 4, 5
        Delta = 1e-6:        Evals needed: 3, 4, 4, 5, 5
        Delta = 1e-6:        Evals needed: 3, 4, 5, 5, 5

    This is as expected, since decreasing delta means Newton's method has
    to approximate the root to a higher degree of accuracy, which generally
    means more iterations needed. 

3) Behavior of x0 in range [1.61, 1.62]

    The actual roots on both sides of 1.62 are about 1.5625 and 1.6667,
    with the midpoint of those roots being 1.6146. Since x in this range is
    close to the midpoint and the derivative at this point, f'(1.62), is
    about -0.86, the first iteration of Newton's method updates the guess
    root, xc, to closer to the first root of the polynomial, which is x =
    1.4706. The method then converges on this root, even though the initial
    guess was closer to the two other roots. 

    However, this is not the case for x0 = 1.61. The derivative at this
    point is about -3.24, which means Newton's method will converge on the
    third root, since the first iteration puts xc close to 1.67, the third
    root.


comment end
%} 


