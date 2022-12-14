

problem = input("Enter Problem Number:\n");

switch problem
    case 1
        %% Homework 6 Problem 1 - The Game of Life
        % Aaron Lee
        % UID: 505540473
        % Objective: Use a Monte Carlo simulation to model John Conway's
        % Game of Life on a 2D grid
        
        % seed the random number generator
        
        rng('shuffle');
        clc; close all; % clear; 
        
        % initialize variables and constants
        
        rows = 150;
        cols = 200;
        trials = 300;
        alive = zeros(trials);
        
        %disp(sum(sum(alive)));
        
        grid = zeros(rows, cols);
        grid_new = zeros(rows, cols);
        
        % fill grid with initial conditions 10% chance of cell starting
        % alive and draw
        
        for i = 1:rows
            for j = 1:cols
                r = rand;
                if r <= 0.1
                    grid(i, j) = 1;
                end
            end
        end
        
        title("Initial Conditions");
        imagesc(grid);
        title("Initial Conditions");
        figure;
        
        % run simulation
        
        for trial = 1:trials
            for i = 1:rows
                for j = 1:cols

                    % define cardinal variables

                    N = i - 1; % north
                    if i == 1
                        N = rows;
                    end

                    S = i + 1; % south
                    if i == rows
                        S = 1;
                    end

                    E = j + 1; % east
                    if j == cols
                        E = 1;
                    end

                    W = j - 1; % west
                    if j == 1
                        W = cols;
                    end

                    % update this generation

                    sum_n = grid(N, W) + grid(N, j) + grid(N, E) + grid(i, W) + grid(i, E) + grid(S, W) + grid(S, j) + grid(S, E);
                    
                    if grid(i, j) == 1
                        if sum_n < 2 || sum_n > 3 % living cell die conditions
                            grid_new(i, j) = 0;
                        end
                    elseif grid(i, j) == 0 && (sum_n == 3) % dead to alive conditions
                        grid_new(i, j) = 1;
                    end
                    
                    % take note of number of alive cells

                    alive(trial) = findSum(grid, rows, cols);
                end
            end
            
            % update old grid and draw this generation
            grid = grid_new;
            title(["Iteration" num2str(trial)]);
            imagesc(grid);
            title(["Iteration" num2str(trial)]);
            drawnow;
        end
        
        % plot cell alive count over time
        figure
        time = 1:1:trials;
        plot(time, alive);
        xlabel("Number of Generations");
        ylabel("Cells Alive");
        title("Number of Cells Alive Over Game of Life Generations");
        
    case 2
        %% Homework 6 Problem 2 - Euler-Bernoulli Beam Bending
        % Aaron Lee
        % UID: 505540473
        % Objective: Model the bending of a beam throughout its length due to an applied force
        % downwards on a point on the beam

        rng('shuffle');
        clc; close all; % clear; 
        
        % initialize constants
        
        l = 1;
        R = 0.013;
        r = 0.011;
        p = -2000;
        d = 0.75;
        E = 70e9;
        I = pi/4 * (R^4 - r^4);
        
        nodes = 20;
        dx = (1/(nodes - 1));
        
        ticks = 0:dx:l;
        A = zeros(nodes, nodes);
        b = zeros(nodes);
        max_disp = [];
        
        % set boundary conditions
        
        A(1, 1) = 1;
        A(nodes, nodes) = 1;
        
        % form constant matrices
        
        % test out different locations of forces (comment out later)
        
        %for d = 0:0.01:l
            for row = 2:nodes - 1
                % fill matrix A with coefficients

                A(row, row - 1) = 1;
                A(row, row) = -2;
                A(row, row + 1) = 1;


                % form right hand side vector of Ay = b, the b vector

                b(row) = dx^2 * M(ticks(row), d, l, p) / (E * I);

            end
            
            % solution to system (displacement vector)
        
            y = A\b;
            results = y(:, 1);
            
            % append location of max displacement to array for optimization
            % problem
            
            [y_max, index] = max(abs(results));
            max_disp(end + 1) = ticks(index);
            
        %end
        
        % find bounds of max displacement for all locations of applied
        % force
        
        max_dispMod = max_disp(max_disp~=0);
        %fprintf("Bounds of Locations of Max Displacement: [%f, %f] m\n", min(max_dispMod), max(max_dispMod));
        
        % plot results of beam bending along rod
        
        plot(ticks, results, 'o-');
        title("Euler-Bernoulli Beam Bending Along Rod");
        xlabel("Position (m)");
        ylabel("Beam Displacement (m)");
        
        % calculate theoritical solution and compare it to graph
        
        fprintf("Max Displacement from Finite Method: %d m\n", -1 * y_max);
        fprintf("Location of Max Displacement Along Rod: %f m\n", max_disp(end));
        
        c = min(d, l - d);
        y_max_t = p * c * (l^2 - c^2)^1.5 / (9 * sqrt(3) * E * I * l);
        
        fprintf("Max Theoretical Displacement: %d m\n", y_max_t);
        fprintf("Error in Displacement: %d m\n", abs(-1 * y_max_t - y_max));
        
    otherwise
        error("Problem Number Invalid");
end

function sum = findSum(grid, rows, cols)
% where grid is a 2D vector, findSum finds the sum of it

    sum = 0;
    for i = 1:rows
        for j = 1:cols
            sum = sum + grid(i, j);
        end
    end
    
end

function m = M(x, d, L, p)
% bending moment of bar function, returns moment as a function of x,
% position along the bar

    if x >= 0 && x <= d
        m = -p * (L - d) * x / L;
    elseif x >= d && x <= L
        m = -p * d *(L - x) / L;
    else
        m = -1;
    end
end