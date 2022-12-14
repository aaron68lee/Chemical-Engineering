

problem = input("Enter Problem Number:\n");

switch problem
    case 1
        %% Homework 1 Problem 1 - Oblate Spheroid Calculations
        % Aaron Lee
        % UID: 505540473
        % Objective: Calculate the surface area of an oblate spheroid given its
        % major and minor radii
        
        % prompt user input
        
        r1 = input("Please input equatorial radius:\n");
        r2 = input("Please input polar radius:\n");
        
        % check validity of user input
        
        if (~(r2 < r1))
            error("Polar radius must be less than equatorial radius");
        end
        
        % calculate SA with equations
        
        y = acos(r2/r1);
        SA = 2*pi * (r1^2 + r2^2/sin(y) * log(cos(y)/(1 - sin(y))));
        SA_approx = 4*pi * ((r1 + r2)/2)^2;
        
        % display SA
        
        fprintf("Surface Area: %9.1f km^2\n", SA);
        fprintf("Approximate SA: %9.1f km^2\n", SA_approx);
        
    case 2
        %% Homework 1 Problem 2 - Nearest Neighbors %%%%%%%%%%%%%%%%%%%%%%
        % Objective: Find the bordering neighbors of a target cell in an N by M
        % grid 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Get user input for dimensions of grid

        n = input("Enter number of rows:\n");
        m = input("Enter number of cols:\n");

        % Error check the validity of dimensions

        mustBeInteger(n) % or use abs(n - floor(n)) < 1e-10
        mustBeInteger(m)

        if (n < 2)
            error("Rows must be at least 2.");
        end
        
        if (m < 2)
            error("Columns must be at least 2.");
        end

        %grid = zeros(n, m);

        % create natural numbers 2D array for debugging purposes
        %{
        for row = 1:n
            for col = 1:m
                grid(row, col) = (col - 1) * n + (row);
            end
        end
        %}

        % Get user input for target cell

        p = input("Enter target cell to check for neighbors:\n");

        mustBeInteger(p);

        % Check target cell validity

        if (p < 1 || p > (n * m))
            error("Cell index out of bounds");
        end

        pRow = mod(p, n);

        if (pRow == 0) % special case for last row
            pRow = 4;
        end

        pCol = ceil(p/n);

        %fprintf("Target cell Row: %i\n", mod(p, n));
        %fprintf("Target cell Col: %i\n", ceil(p/n));

        % create empty array to add neighbors to

        neighbors = [];

        % check back 1 column

        if(pCol - 1 > 0) 
            if (pRow - 1 > 0)
                neighbors(end + 1) = p - n - 1;
            end
            neighbors(end + 1) = p - n;
            if (pRow + 1 <= n)
                neighbors(end + 1) = p - n + 1;
            end
        end

        % check current column

        if (pRow - 1 > 0)
            neighbors(end + 1) = (p - 1);
        end

        if (pRow + 1 <= n)
            neighbors(end + 1) = (p + 1);
        end

        % check forward 1 column

        if (pCol + 1 <= m) 
            if (pRow - 1 > 0)
                neighbors(end + 1) = (p + n - 1);
            end

            neighbors(end + 1) = (p + n);

            if (pRow + 1 <= n)
                neighbors(end + 1) = (p + n + 1);
            end
        end

        % print all neighbors

        fprintf('For grid with %d rows and %d columns: \n', n, m)
        
        a = 'Cell ID:   ';
        b = num2str(p);

        disp([a b]);

        fprintf("Neighbors: ");

        for e = 1:length(neighbors)
            fprintf("%i ", neighbors(e));
            %neighbor_list = [neighbor_list num2str(neighbors(e))];
        end
        
        fprintf("\n");
        %{
        a = 'Neighbors: ';
        b = neighbor_list;
        disp([a b]);
        %}
    otherwise
        error("Problem Number Invalid");
end
