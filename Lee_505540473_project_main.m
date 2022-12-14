%% Final Project - Spinodal Decomposition
% Aaron Lee
% UID: 505540473
% Objective: Use phase field theory to analyze the spinodal decomposition
% evolution of various alloys over time
        
clc; clear; close all;

% initialize constants 

a = 0.1;
b = 1;
gamma = 1;
D = 3;
h = 1;

dt = [1e-3, 1e-4, 1e-4];
t0 = 0;
tf = 20;
rows = 150;
cols = 100;
fps = 30;

steps = ceil(tf./dt);

% 9 point stencil first order Runge-Kutta
% 9 point stencil fourth order Runge-Kutta

% reinitialize grid
% generate random grid

grid = genGrid(rows, cols);
old_grid = grid;

title("Initial Conditions");
xlim([1 100])
ylim([1 150])
imagesc(old_grid);
colorbar;
title("Initial Conditions");

figure;
        

%% Run simulation of spinodal decomposition using Runge-Kutta and stencil approximations

for trial = 1:3
    % initialize stencil variables per trial
    if trial == 1
        stencil = 5;
    else
        stencil = 9;
    end
    
    if trial == 3
        suffix = "_4th_Order.mp4";
    else
        suffix = "_1st_Order.mp4";
    end
    
    % initialize video file
    
    name = num2str(stencil) + "_Stencil" + suffix;
    video_file = VideoWriter(name, 'MPEG-4');
    video_file.FrameRate = fps;
    open(video_file);

    % reset to original distribution for comparisons between stencils
    
    grid = old_grid;
    
    for k = 1:steps(trial)
        
        % update this generation

        % current time is k*steps(trial)

        if trial ~= 3 % first order runge-kutta
            
            L = Laplacian_2D(grid, h, stencil);
            grid_temp = b^4 * grid.^3 - a*(b^2)*grid - gamma*L;

            c1 = dt(trial) * D * Laplacian_2D(grid_temp, h, stencil);
            grid_new = grid + c1;
            
        else % fourth order runge-kutta
            
            % do fourth order approximation here
            
            L = Laplacian_2D(grid, h, stencil);
            grid_temp = b^4 * grid.^3 - a*(b^2)*grid - gamma*L;

            c1 = dt(trial) * D * Laplacian_2D(grid_temp, h, stencil);
            c2 = dt(trial) * D * Laplacian_2D(grid + c1/2, h, stencil);
            c3 = dt(trial) * D * Laplacian_2D(grid + c2/2, h, stencil);
            c4 = dt(trial) * D * Laplacian_2D(grid + c3, h, stencil);
            grid_new = grid + c1/6 + c2/3 + c3/3 + c4/6;
            
        end
              
        % update old grid and draw this generation with 30 fps
        
        if (trial == 1 && mod(k * dt(trial), 1/30) <= 7e-4) || ((trial == 2 || trial == 3) && mod(k * dt(trial), 1/30) <= 7e-5) % produce 30 fps framerate
            title(["Trial" num2str(trial) " Iteration" num2str(k)]);
            imagesc(grid_new);
            colorbar;
            title(["Trial" num2str(trial) " Iteration" num2str(k)]);
            drawnow;
            
            % save frame to video
            
            frame = getframe(gcf);
            writeVideo(video_file, frame);
            
        end
        
        % print plot every 5 seconds
        
        if mod(k * dt(trial), 5) <= 1e-5 
            title(["Trial" num2str(trial) " Timestamp(s): " num2str(k * dt(trial))]);
            imagesc(grid_new);
            colorbar;
            title(["Trial" num2str(trial) " Timestamp(s): " num2str(k * dt(trial))]);
            drawnow;
            figure;
        end
        
        grid = grid_new;
        
    end
    close(video_file);
end
    

%% define Lapalcian function

function v = Laplacian_2D(u, h, stencil)
% performs a laplacian on input matrix, u, with discretized grid spacing of
% h using stencil, which is either 5 or 9-point approximation
% returns v, the solution of the laplacian

    % define constants and initialize solution matrix
    
    rows = 150;
    cols = 100;
    
    v = zeros(rows, cols);
    
    %{
    switch stencil
        case 5
        case 9
        otherwise
            error("Stencil Number Invalid");
    end
    %}
    
    for i = 1:rows
        for j = 1:cols

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
            
            if stencil == 5
                v(i, j) = 1/(h^2) * ( u(i, W) + u(i, E) + u(N, j) + u(S, j) - 4*u(i, j) );
            elseif stencil == 9
                v(i, j) = 1/(6 * h^2) * ( u(N, W) + u(N, E) + u(S, W) + u(S, E) + 4*u(i, W) + 4*u(i, E) + 4*u(N, j) + 4*u(S, j) - 20*u(i, j) );
            else
                error("Stencil Number Invalid");
            end
            
        end
    end
            
end

%% define Grid Generation function

function grid = genGrid(r, c)
% generates a random grid with different material phases with dimensions of
% r rows and c cols
    
    grid = ones(r, c);

    for i = 1:r
        for j = 1:c
            r = rand;
            if r <= 0.5
                grid(i, j) = -1;
            end
        end
    end
end