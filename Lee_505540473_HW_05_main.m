

problem = input("Enter Problem Number:\n");

switch problem
    case 1
        %% Homework 5 Problem 1 - The Shared Birthday Problem
        % Aaron Lee
        % UID: 505540473
        % Objective: The Shared Birthday problem aims to use a Monte Carlo 
        % simulation to determine the number of people required in a group 
        % to have a pair with birthdays less than seven days apart. 
        
        rng('shuffle');
        clc; close all; clear; 
        
        trials = 1e4;
        results = zeros(trials, 1);
        group = [];
        count = 0;
        
        for i = 1:trials
            count = 0;
            group = [];
            while ~birthdayTest(group) % while no two people have birthdays in same week
                group(end + 1) = ceil(365 * rand());
                count = count + 1;
            end
            results(i) = count;
        end
        
        fprintf("Median Number of People = %2i\n", median(results));
        histogram(results);
        title("Distribution of Number of People Needed for Shared Weekly Birthdays");
        xlabel('Birthday Day');
        ylabel('Frequency');
        
    case 2
        %% Homework 5 Problem 2 - Random Walk Collisions
        % Aaron Lee
        % UID: 505540473
        % Objective: 

        rng('shuffle');
        clc; close all; % clear; 
        
        % initialize constants and initial conditions
        
        p = 0.2;
        rows = 11; cols = 11;
     
        grid = zeros(rows, cols);
        
        %{
        posA = [ceil(rows/2) 1];
        posB = [ceil(rows/2) cols];
        positionsA = [];
        positionsB = [];
        
        r1 = posA(1);
        c1 = posA(2);
        r2 = posB(1);
        c2 = posB(2);
        %}
        
        moves = 0;
        maxMoves = 1000;
        
        
        %%%%%%%%%%%%%%%%%%% single walker simulation %%%%%%%%%%%%%%%%%%
        
        %{
        while moves < maxMoves % still within grid boundaries and less than max moves
            
            % log the current position
            
            positionsA = [positionsA; posA];
            
            % Move particla A
            
            r = rand;

            if r <= p
                % Move North
                if(posA(1) >= 2)
                    r1 = r1 - 1;
                    % Update Position
                    posA = [r1, c1];
                end
            elseif p < r && r <= 2*p
                % Move East
                if(posA(2) < cols)
                    c1 = c1 + 1;
                    % Update Position
                    posA = [r1, c1];
                end
            elseif 2*p < r && r <= 3*p
                % Move South
                if(posA(1) < rows)
                    r1 = r1 + 1;
                    % Update Position
                    posA = [r1, c1];
                end
            elseif 3*p < r && r <= 4*p
                % Move West
                if(posA(2) >= 2)
                    c1 = c1 - 1;
                    % Update Position
                    posA = [r1, c1];
                end
            elseif (1 - 4*p) < r
                % Stay put
            end

            moves = moves + 1;
            
        end
        %}
        
        %%%%%%%%%%%%%%%%%% two-walker collision simulation %%%%%%%%%%%%%%%%%%
        
        % initialize initial conditions
        
        collision = false;
        trials = 5000;
        results = zeros(trials, 1);
        
        % plot initial conditions
        
        posA = [6 3];
        posB = [6 9];
            
        grid(posA(1), posA(2)) = 0.5;
        grid(posB(1), posB(2)) = 1;
        imagesc(grid);
        figure;
        
        
        for i = 1:trials
            
            % reset positions and counter variables
            
            grid = zeros(rows, cols);
            
            moves = 0;
            posA = [6 1];
            posB = [6 11];
            positionsA = [];
            positionsB = [];

            r1 = posA(1);
            c1 = posA(2);
            r2 = posB(1);
            c2 = posB(2);
            
            % draw initial conditions of both particles
            
            grid(posA(1), posA(2)) = 0.5;
            grid(posB(1), posB(2)) = 1;
            
            collision = false;
            
            while ~collision && moves < maxMoves % no collision exists and move still within max move range
                
                % log the current positions for debugging
                %{
                positionsA = [positionsA; posA];
                positionsB = [positionsB; posB];
                %}
                
                % Move particla A

                r = rand;

                if r <= p
                    % Move North
                    if(posA(1) >= 2)
                        r1 = r1 - 1;
                        % Update Position
                        grid(posA(1), posA(2)) = 0;
                        posA = [r1, c1];
                        grid(posA(1), posA(2)) = 0.5;
                    end
                elseif p < r && r <= 2*p
                    % Move East
                    if(posA(2) < cols)
                        c1 = c1 + 1;
                        % Update Position
                        grid(posA(1), posA(2)) = 0;
                        posA = [r1, c1];
                        grid(posA(1), posA(2)) = 0.5;
                    end
                elseif 2*p < r && r <= 3*p
                    % Move South
                    if(posA(1) < rows)
                        r1 = r1 + 1;
                        % Update Position
                        grid(posA(1), posA(2)) = 0;
                        posA = [r1, c1];
                        grid(posA(1), posA(2)) = 0.5;
                    end
                elseif 3*p < r && r <= 4*p
                    % Move West
                    if(posA(2) >= 2)
                        c1 = c1 - 1;
                        % Update Position
                        grid(posA(1), posA(2)) = 0;
                        posA = [r1, c1];
                        grid(posA(1), posA(2)) = 0.5;
                    end
                elseif 4*p < r
                    % Stay put
                    %fprintf("%f\n", 1-4*p);
                end

                % Move particle B
                % comment this part to fix particle B position
                
                
                r = rand;

                if r <= p
                    % Move North
                    if(posB(1) >= 2)
                        r2 = r2 - 1;
                        % Update Position
                        grid(posB(1), posB(2)) = 0;
                        posB = [r2, c2];
                        grid(posB(1), posB(2)) = 1;
                    end
                elseif p < r && r <= 2*p
                    % Move East
                    if(posB(2) < cols)
                        c2 = c2 + 1;
                        % Update Position
                        grid(posB(1), posB(2)) = 0;
                        posB = [r2, c2];
                        grid(posB(1), posB(2)) = 1;
                    end
                elseif 2*p < r && r <= 3*p
                    % Move South
                    if(posB(1) < rows)
                        r2 = r2 + 1;
                        % Update Position
                        grid(posB(1), posB(2)) = 0;
                        posB = [r2, c2];
                        grid(posB(1), posB(2)) = 1;
                    end
                elseif 3*p < r && r <= 4*p
                    % Move West
                    if(posB(2) >= 2)
                        c2 = c2 - 1;
                        % Update Position
                        grid(posB(1), posB(2)) = 0;
                        posB = [r2, c2];
                        grid(posB(1), posB(2)) = 1;
                    end
                elseif 4*p < r
                    % Stay put
                end
                
                
                
                % Check for a collision
                if posA == posB
                    collision = true;
                end

                moves = moves + 1;
                
                % draw the visual of particles moving
                
                imagesc(grid);
                drawnow;
                
            end
            results(i) = moves;
            
        end
        
        fprintf("Median = %2i\n", median(results));
        %histogram(results);
        
    otherwise
        error("Problem Number Invalid");
end