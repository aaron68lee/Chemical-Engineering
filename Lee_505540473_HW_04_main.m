

problem = input("Enter Problem Number:\n");

switch problem
    case 1
        %% Homework 4 Problem 1 - The Split and Average Problem
        % Aaron Lee
        % UID: 505540473
        % Objective: Create a smoothened curve or surface from an input set
        % of points by "splitting" -- creating intermediary points-- and
        % "averaging" them with weights for a smoothened effect.
        
        % initialize weights and point arrays
       
        x = [0, 0, 1, 1];
        y = [0, 1, 0, 1];
        
        x = [0, 0.2, 0.4, 0.7, 1, 0.3, 0.8];
        y = [0, 0.5, 0.8, 0.9, 1, 0.2, 0.3];
        
        w = [100, 1, -50];
        
        % split and average recursively
        
        dx = 1;
        dy = 1;
        
        % plot initial set of points
        
        plot(x, y, 'bo');
        figure
        count = 0; % count average number of iterations until convergence
        
        while (max(sqrt(dx.^2 + dy.^2)) >= 1e-3 && count <= 15)
            
            % create split array of points
            
            xs = splitPts(x);
            ys = splitPts(y);
            
            % weight average splitted arrays
            
            xa = averagePts(xs, w);
            ya = averagePts(ys, w);
            x = xa;
            y = ya;
            
            % determine max displacement for convergence test
            
            dx = xa - xs;
            dy = ya - ys;
            count = count + 1;
        end
        
        % plot final set of points
        
        plot(x, y, 'bo');
        fprintf("Iterations Until Convergence: %i\n", count);
        fprintf("Number of Points: %i\n", length(xa));
        
    case 2
        %% Homework 4 Problem 2 - Runge-Kutta Radioactivity
        % Aaron Lee
        % UID: 505540473
        % Objective: The goal of the Runge-Kutta method is to approximate a
        % function using different degrees of accuracies and orders, 
        % specifically, exponential radioactive decay, and
        % comparing it to the exact differential solution
        
        clc; clear;
        
        dt = [1, 0.1, 0.01];
        t0 = 0;
        tf = 15;
        y = 1;
        hl = 2.45;
        steps = ceil(tf./dt);
        
        t1 = linspace(t0, tf, steps(1));
        t2 = linspace(t0, tf, steps(2));
        t3 = linspace(t0, tf, steps(3));
        
        y1 = [];
        y2 = [];
        y4 = [];
        
        fprintf("   dt       RK1       RK2       RK4\n");
        
        % plot exact solution
        
        hold on
        
        % first order approximation time step 1
        
        y = 1;
        y_sol1 = y * exp((-log(2)/hl) * t1);
        
        for k = 1:steps(1)
           y1(end + 1) = y; 
           y = advanceRK(y, dt(1), 1);
        end
        
        title("Carbon-15 Decay Using Timestep 1.00s");
        plot(t1, y_sol1, 'DisplayName', 'Exact Solution');
        
        xlabel("Time (s)");
        ylabel("Amount of Carbon-15 (arbitrary units)");
        legend
        
        plot(t1, y1, 'DisplayName', "First Order Runge-Kutta Approximation");
        
        % second order approximation time step 1
        
        y = 1;
        for k = 1:steps(1)
           y2(end + 1) = y;
           y = advanceRK(y, dt(1), 2); 
        end
        
        plot(t1, y2, 'DisplayName', "Second Order Runge-Kutta Approximation");
        
        % fourth order approximation time step 1
        
        y = 1;
        for k = 1:steps(1)
           y4(end + 1) = y;
           y = advanceRK(y, dt(1), 4); 
        end
        
        plot(t1, y4, 'DisplayName', "Fourth Order Runge-Kutta Approximation");
        
        % find and print time step 1 errors
        
        err11 = mean(abs(y1 - y_sol1));
        err12 = mean(abs(y2 - y_sol1));
        err14 = mean(abs(y4 - y_sol1));
        
        fprintf("%1.2f:  %1.2e  %1.2e  %1.2e\n", dt(1), err11, err12, err14);
        
        hold off
        figure
        
        % time step 2
        
        y = 1;
        y1 = [];
        y2 = [];
        y4 = [];
        
        hold on
        
        y_sol2 = y * exp(-log(2)/hl * t2);
        title("Carbon-15 Decay Using Timestep 0.10s");
        plot(t2, y_sol2, 'DisplayName', 'Exact Solution');
        
        xlabel("Time (s)");
        ylabel("Amount of Carbon-15 (arbitrary units)");
        legend
        
        % first order approximation time step 2
        
        y = 1;
        for k = 1:steps(2)
           y1(end + 1) = y; 
           y = advanceRK(y, dt(2), 1);
        end
        
        plot(t2, y1, 'DisplayName', "First Order Runge-Kutta Approximation");

        % second order approximation time step 2
        
        y = 1;
        for k = 1:steps(2)
           y2(end + 1) = y;
           y = advanceRK(y, dt(2), 2); 
        end
        
        plot(t2, y2, 'DisplayName', "Second Order Runge-Kutta Approximation");
        
        % fourth order approximation time step 2
        
        y = 1;
        for k = 1:steps(2)
           y4(end + 1) = y;
           y = advanceRK(y, dt(2), 4); 
        end
        
        plot(t2, y4, 'DisplayName', "Fourth Order Runge-Kutta Approximation");
        
        % find and print time step 2 errors
        
        err21 = mean(abs(y1 - y_sol2));
        err22 = mean(abs(y2 - y_sol2));
        err24 = mean(abs(y4 - y_sol2));
        
        fprintf("%1.2f:  %1.2e  %1.2e  %1.2e\n", dt(2), err21, err22, err24);
        
        hold off
        figure
        
        % time step 3
        
        y1 = [];
        y2 = [];
        y4 = [];
        y = 1;
        
        hold on
        
        y_sol3 = y * exp(-log(2)/hl * t3);
        title("Carbon-15 Decay Using Timestep 0.01s");
        plot(t3, y_sol3, 'DisplayName', 'Exact Solution');
        
        xlabel("Time (s)");
        ylabel("Amount of Carbon-15 (arbitrary units)");
        legend
        
        % first order approximation time step 3
        t_y = [];
        y = 1;
        for k = 1:steps(3)
           y1(end + 1) = y; 
           y = advanceRK(y, dt(3), 1);
        end
        
        plot(t3, y1, 'DisplayName', "First Order Runge-Kutta Approximation");

        % second order approximation time step 3
        
        y = 1;
        for k = 1:steps(3)
           y2(end + 1) = y;
           y = advanceRK(y, dt(3), 2); 
        end
        
        plot(t3, y2, 'DisplayName', "Second Order Runge-Kutta Approximation");
        
        % fourth order approximation time step 3
        
        y = 1;
        for k = 1:steps(3)
           y4(end + 1) = y;
           y = advanceRK(y, dt(3), 4); 
           t_y = [t_y y];
        end
        
        plot(t3, y4, 'DisplayName', "Fourth Order Runge-Kutta Approximation");
        
        % find and print time step 3 errors
        
        err31 = mean(abs(y1 - y_sol3));
        err32 = mean(abs(y2 - y_sol3));
        err34 = mean(abs(y4 - y_sol3));
        
        fprintf("%1.2f:  %1.2e  %1.2e  %1.2e\n", dt(3), err31, err32, err34);
        fprintf("Percent of Carbon-15 left after %i seconds: %2.2f percent\n", tf, y4(end) * 100); % percent left as scientific notation
        
        hold off
        
    otherwise
        error("Problem Number Invalid");
end

%%%%%%%%%%%%%%%%%%%%%%% Averaging function %%%%%%%%%%%%%%%%%%%%%%%%%

function [ xa ] = averagePts( xs, w )
    % the averagePts function takes an input vector xs with the points to
    % average and an input vector w, containing 3 weights
    % averagePts returns a vector, xa, the same size as xs, containing
    % weighted-averages of each group of 3 centered on each element of xs.
    if sum(w) == 0
        error("Sum of Weights Cannot be 0");
    else
        nw = w/sum(w);
    end
    
    size = length(xs);
    xa = (size);

    for i = 1:size
        if i == 1
            xa(i) = nw(1)*xs(size) + nw(2)*xs(i) + nw(3)*xs(i+1);
        elseif i == size
            xa(i) = nw(1)*xs(i-1) + nw(2)*xs(i) + nw(3)*xs(1);
        else
            xa(i) = nw(1)*xs(i-1) + nw(2)*xs(i) + nw(3)*xs(i+1);
        end
    end
    
    % print return vector for debugging 
    %{
    for k = 1:size
            fprintf("%.4f, ", xa(k));
    end
    %}  
    
end

%%%%%%%%%%%%%%%%%%%%%%% Splitting function %%%%%%%%%%%%%%%%%%%%%%%%%

function [ xs ] = splitPts(x)
        % this function returns a vector containing the original points in
        % vector, x, with midpoints inserted between elements of x
        
        size = length(x);
        xs = (size * 2);
        
        for i = 1:size
            % fill odd positions with original elements
            xs(i * 2 - 1) = x(i);
            
            % special case for last element
            if i ~= size
                xs(i * 2) = (x(i) + x(i + 1))/2;
            else
                xs(i * 2) = (x(i) + x(1))/2;
            end
            
        end
        
        % print return vector for debugging 
        %{
        for k = 1:size * 2
            fprintf("%.4f, ", xs(k));
        end
        %}
end
        
%%%%%%%%%%%%%%%%%%%%%%% runge-kutta function %%%%%%%%%%%%%%%%%%%%%%%%%

function [ y ] = advanceRK(y, dt, method)
% this function advances a discretized solution by one time step
% where y is the current amount of carbon-15
% dt is the time-step size
% and method is either 1, 2, or 4, determining the order of approximation
    hl = 2.45;
    
    switch method

        case 1
            c1 = dt * y * -log(2)/hl;
            y_new = y + c1;
            y = y_new;
        case 2
            c1 = dt * y * -log(2)/hl;
            c2 = dt * (y + c1/2) * -log(2)/hl;
            y_new = y + c2;
            y = y_new;
        case 4
            c1 = dt * y * -log(2)/hl;
            c2 = dt * (y + c1/2) * -log(2)/hl;
            c3 = dt * (y + c2/2) * -log(2)/hl;
            c4 = dt * (y + c3) * -log(2)/hl;
            y_new = y + c1/6 + c2/3 + c3/3 + c4/6;
            y = y_new;
        otherwise
            error("Approximation Order Invalid.");
    end

end
