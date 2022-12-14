

problem = input("Enter Problem Number:\n");

switch problem
    case 1
        %% Homework 2 Problem 1 - Three Species Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Aaron Lee
        % UID: 505540473
        % Objective: Simulate a three-species predatory-prey model defined
        % by initial population conditions over a period of time modeled by
        % given equations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % set initial conditions 
        
        % x_new denotes the new state, while x denotes the old state
        
        
        x = 2;
        y = 2.49;
        z = 1.5;
        
        dt = 0.005;
        tf = 12;
        steps = ceil(tf/dt);
        t0 = 0;
        test_size = 0.1;
        zero = 0.01;
        allLived = 0;
        
        % print initial conditions
        
        fprintf("Time    X    Y    Z\n");
        fprintf(" %2.1f %1.2f %1.2f %1.2f\n", t0, x, y, z);
        
        % time how long it takes to simulate using the model 
        
        %{
        tic 
        for a = 1:test_size:30
            for b = 1:test_size:30
                for c = 1:test_size:30
                    
                    a0 = a;
                    b0 = b;
                    c0 = c;
                    
                    x = a;
                    y = b;
                    z = c;
                    
                    for k = 1:steps
            
                        % propogate new population values

                        x_new = dt * (0.75 * x * (1-x/20) - 1.5*x*y - 0.5*x*z) + x;
                        y_new = dt * (y * (1-y/25) - 0.75*x*y - 1.25*y*z) + y;
                        z_new = dt * (1.5 * z * (1-z/30) - x*z - y*z) + z;

                        % update population values with new values

                        x = x_new;
                        y = y_new;
                        z = z_new;

                    end
                    
                    if(x >= zero && y >= zero && z >= zero)
                        fprintf("Final pops: %1.2f %1.2f %1.2f\n", x, y, z);
                        fprintf("      Initial Conditions: %1.2f %1.2f %1.2f\n", a0, b0, c0);
                        allLived = allLived + 1;
                    end
                    
                end
            end
        end
        
        fprintf("Number Cases Where All Lived: %i\n", allLived);
        %}
        
        
        for k = 1:steps
            
            % propogate new population values
            
            x_new = dt * (0.75 * x * (1-x/20) - 1.5*x*y - 0.5*x*z) + x;
            y_new = dt * (y * (1-y/25) - 0.75*x*y - 1.25*y*z) + y;
            z_new = dt * (1.5 * z * (1-z/30) - x*z - y*z) + z;
            
            % print popoulation data every 0.5 seconds
            
            if(mod(dt * k, 0.5) < 1e-10)
                %fprintf("%2.1f %1.2f %1.2f %1.2f\n", 0.005*k, x, y, z);
                if(dt*k < 10)
                    fprintf(" %2.1f %1.2f %1.2f %1.2f\n", dt*k, x_new, y_new, z_new);
                else
                    fprintf("%2.1f %1.2f %1.2f %1.2f\n", dt*k, x_new, y_new, z_new);
                end
            end
            
            % update population values with new values
            
            x = x_new;
            y = y_new;
            z = z_new;
            
        end
        %
        toc
        
    case 2
        %% Homework 2 Problem 2 - The Pocket Change Problem
        % Objective: Find the number of coins needed for a given amount of
        % change as well as average number of coins needed for all change
        % values for different denominations
        
        % initialize coin values
        
        totalcount = 0;
        QUARTER = 25;
        DIME = 10;
        NICKEL = 5;
        PENNY = 1;
        
        interval = 1;
        trials = 0;
        
        % loop through all possible change values
        
        for i = 0:interval:99
            
            change = i;
            count = 0;
            
            % find number of each coin needed for this iteration of change
            % value
            
            while change > 1e-10
                if(change >= QUARTER)
                    change = change - 25;
                    count = count + 1;
                elseif(change >= DIME)
                    change = change - 10;
                    count = count + 1;
                elseif(change >= NICKEL)
                    change = change - 5;
                    count = count + 1;
                else
                    change = change - 1;
                    count = count + 1;
                end
            end
            
            % for debugging only
            
            %fprintf("Change: %i, Coins: %i\n", i, count);
            
            % add to totals
            
            trials = trials + 1;
            totalcount = totalcount + count;
        end
        
        % print results
        
        fprintf("\nAverage Number of Coins = %.2f\n", totalcount/trials);
        %fprintf("# Trials: %i\n", trials);
        
    otherwise
        error("Problem Number Invalid");
end
