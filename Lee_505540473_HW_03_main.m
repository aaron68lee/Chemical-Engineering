

problem = input("Enter Problem Number:\n");

switch problem
    case 1
        %% Homework 3 Problem 1 - The Pendulum Physics Problem
        % Aaron Lee
        % UID: 505540473
        % Objective: 
        
        % set initial conditions 
        
        % x_new denotes the new state, while x denotes the old state
        
        L = 1;
        theda = pi/3;
        t0 = 0;
        w = 0;
        g = 9.81;
        
        position = [];
        velocity = [];
        accel = [];
        energy = [];
        
        position2 = [];
        velocity2 = [];
        accel2 = [];
        energy2 = [];
        
        dt = 0.0005;
        tf = 20;
        steps = ceil(tf/dt);
        
        % print initial conditions
        %{
        fprintf("Time    Position    Angular Velocity    Acceleration    Energy\n");
        fprintf(" %2.1f       %.2f                %.2f             %.2f       %.2f\n", t0, theda, 0, -8.5, nrg(theda0, 0, L));
        %}
        % loop using explicit euler method
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for k = 1:steps
            
            % propogate new population values
            
            theda_new = dt * w + theda;
            w_new = dt * (-g/L * sin(theda)) + w;
            
            % add values of position and velocity to arrays
            
            position(end + 1) = theda;
            velocity(end + 1) = w;
            accel(end + 1) = -g/L * sin(theda);
            energy(end + 1) = nrg(theda, w, L);
            
            % print position, velocity, and energy data every 0.5 seconds
            %{
            if(mod(dt * k, 0.5) < 1e-10)
                if(dt*k < 10)
                    fprintf(" %2.1f       %.2f                %.2f             %.2f       %.2f\n", dt*k, theda_new, w_new, -g/L*sin(theda), nrg(theda_new, w_new, L));
                else
                    fprintf("%2.1f       %.2f                %.2f             %.2f       %.2f\n", dt*k, theda_new, w_new, -g/L*sin(theda), nrg(theda_new, w_new, L));
                end
            end
            %}
            % update old values with new values
            
            theda = theda_new;
            w = w_new;
           
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        L = 1;
        theda = pi/3;
        t0 = 0;
        w = 0;
        
        %{
        fprintf("\n\nTime    Position    Angular Velocity    Acceleration    Energy\n");
        fprintf(" %2.1f       %.2f                %.2f             %.2f       %.2f\n", t0, theda, 0, 0, nrg(theda0, 0, L));
        %}
        
        % loop using semi-implicit euler method
        
        for k = 1:steps
            
            % propogate new population values
            
            w_new = dt * (-g/L * sin(theda)) + w;
            theda_new = dt * w_new + theda;
            
            % add values of position and velocity to arrays
            
            position2(end + 1) = theda;
            velocity2(end + 1) = w;
            accel2(end + 1) = -g/L * sin(theda);
            energy2(end + 1) = nrg(theda_new, w_new, L);
            
            % print position, velocity, and energy data every 0.5 seconds
            %{
            if(mod(dt * k, 0.5) < 1e-10)
                if(dt*k < 10)
                    fprintf(" %2.1f       %.2f                %.2f             %.2f       %.2f\n", dt*k, theda_new, w_new, -g/L*sin(theda), nrg(theda_new, w_new, L));
                else
                    fprintf("%2.1f       %.2f                %.2f             %.2f       %.2f\n", dt*k, theda_new, w_new, -g/L*sin(theda), nrg(theda_new, w_new, L));
                end
            end
            %}
            % update old values with new values
            
            theda = theda_new;
            w = w_new;
           
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        time = linspace(t0, tf, steps);
        
        % plot explicit pendulum graph
        
        hold on;
        plot(time, position, 'DisplayName', 'Angular Position');
        plot(time, velocity, 'DisplayName', 'Angular Velocity');
        plot(time, accel, 'DisplayName', 'Angular Acceleration');
        hold off;
        
        xlabel("Time (s)");
        ylabel("Position (rads)");
        title("Pendulum Swings Over Time | Explicit Euler");
        axis([0 20 -2*pi 2*pi]);
        ylim([min(accel), max(accel)]);
        legend
        
        figure
        
        % plot explicit pendulum energy 
        graph
        plot(time, energy, 'DisplayName', 'Energy');
        
        xlabel("Time (s)");
        ylabel("Energy (J)");
        title("Pendulum Energy Over Time | Explicit Euler");
        axis([0 20 -2*pi 2*pi]);
        ylim([min(energy), max(energy)]);
        legend
        
        figure
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % plot semi-implicit pendulum graphs
        
        hold on;
        plot(time, position2, 'DisplayName', 'Angular Position');
        plot(time, velocity2, 'DisplayName', 'Angular Velocity');
        plot(time, accel2, 'DisplayName', 'Angular Acceleration');
        hold off;
        
        xlabel("Time (s)");
        ylabel("Position (rads)");
        title("Pendulum Swings Over Time | Semi-implicit Euler");
        axis([0 20 -2*pi 2*pi]);
        ylim([min(accel2), max(accel2)]);
        legend
        
        figure
        
        % plot semi-implicit pendulum energy 
        
        graph
        plot(time, energy2, 'DisplayName', 'Energy');
        
        xlabel("Time (s)");
        ylabel("Energy (J)");
        title("Pendulum Energy Over Time | Semi-implicit Euler");
        %axis([0 20 -2*pi 2*pi]);
        %ylim([min(energy2), max(energy2)]);
        legend
        
    case 2
        %% Homework 3 Problem 2 - DNA Analysis
        % Objective: 
        
        % load in dna data and initialize variables
        
        load('chr1_sect.mat');
        BASES = length(dna);
        startP = 0;
        endP = 0;
        proteins = [];
        TAA = 0;
        TAG = 0;
        TGA = 0;
        
        for i = 1:3:BASES - 2
            
            % check for start codons
            
            if startP == 0
                if dna(i) == 1 && dna(i+1) == 4 && dna(i+2) == 3
                    startP = i; % set start to current position
                end
            else % check for stop codons
                if (dna(i) == 4 && dna(i+1) == 1 && dna(i+2) == 1) || (dna(i) == 4 && dna(i+1) == 1 && dna(i+2) == 3) || (dna(i) == 4 && dna(i+1) == 3 && dna(i+2) == 1)
                   
                    % count how many of each type of stop codon encountered
                    
                   if(dna(i) == 4 && dna(i+1) == 1 && dna(i+2) == 1)
                       TAA = TAA + 1;
                   elseif (dna(i) == 4 && dna(i+1) == 1 && dna(i+2) == 3)
                       TAG = TAG + 1;
                   else
                       TGA = TGA + 1;
                   end
                   
                   % append this protein length to protein array
                   
                   proteins(end + 1) = i - startP + 3;
                   startP = 0;
                end
            end
        end
        
        fprintf("Total Protein-Coding Segments: %i\n", length(proteins));
        fprintf("Average Length: %.2f\n", mean(proteins));
        fprintf("Maximum Length: %i\n", max(proteins));
        fprintf("Minimum Length: %i\n\n", min(proteins));
        fprintf("Bases used for proteins: %i\n", sum(proteins));
        fprintf("Percent Proteins: %2.2f\n", sum(proteins)/BASES * 100);
        fprintf("TAA: %i, TAG: %i, TGA: %i\n", TAA, TAG, TGA);
        disp(BASES);
        
    otherwise
        error("Problem Number Invalid");
        % define energy function at given time step
        
        
end

function e = nrg(theda, w, L)
    h = L*(1-cos(theda));
    e = 9.81 * h + 0.5*(L*w)^2;
end