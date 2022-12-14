
%% Extra Credit Problem 1 - Ranked Choice Vote
% Aaron Lee
% UID: 505540473
% Objective: Simulate a ranked choice vote, a voting system where only a
% voter's highest ranked vote is considered and lowest ranking candidates
% are systematically eliminated

clc; clear; close all;

% load data files

%load('votes1.mat');
load('votes2.mat');

%{
votes = [1, 2, 3; 3, 2, 1; 1, 3, 2; 4, 3, 2];
%}
%votes = [3, 1, 2, 4; 2, 4, 3, 1; 2, 1, 3, 4; 4, 1, 3, 2; 2, 1, 3, 4; 1, 3, 4, 2; 2, 1, 3, 4; 4, 1, 3, 2; 3, 1, 4, 2];
 

candidates = length(votes(1, :));
header = 1:1:candidates;
win = false;
winner = 0;
count = 0;

% run voting system

%{
disp("before...");
disp(votes);
%}


% Modify for weighted voting system
%{
totals = zeros(1, candidates);

for i = 1:length(votes)
    for j = 1:candidates
        totals(votes(i, j)) = totals(votes(i, j)) + (candidates - j + 1);
        %votes(:, i) = votes(:, i) * (candidates - i + 1);
    end
end

fprintf("Totals:     " + num2str(totals) + "\n");
[M, indexMax] = max(totals);
fprintf("Winning Candidate: %i\n", indexMax);
%}

% print header

fprintf("                    "); % + num2str(header) + "\n");
for i = 1:candidates
    fprintf("%4.f  ", header(i));
end
fprintf("\n");

while ~win % majority winner not found
    
    % calculate current round totals
    
    favs = votes(:, 1);
    
    for i = 1:candidates
        temp = favs(favs == i);
        totals(i) = length(temp);
    end
    
    % calculate most and least popular 
    
    [maximum, indexMax] = max(totals);
        
    if maximum/sum(totals) > 0.5
        win = true;
        winner = indexMax;
    else
      
        %[min, indexMin] = min(totals);
        
        % find least popular nonzero candidate
        
        min = realmax;
        indexMin = 0;
        for i = 1:candidates
            if totals(i) ~= 0 && totals(i) < min
                min = totals(i);
                indexMin = i;
            end
        end
        
        % remove least popular candidate (that is nonzero)
        
        votes = removeCandidate(votes, indexMin);

    end
    count = count + 1;
    
    % disp current round results 
    
    fprintf("Round %i Totals:     ", count); % + num2str(totals) + "\n", count);
    for i = 1:candidates
        fprintf("%4.f  ", totals(i));
    end
    fprintf("\n");
    %disp(votes);
end

fprintf("Winning Candidate:  %i\n", winner);
    
%% function implementations

function votes = removeCandidate(votes,losingCandidate)
% votes is a 2d array of votes, with each voter's choices representing a single row
% removeCandidate removes all instances of the losingCandidate number from the votes array
    
    rows = length(votes);
    cols = length(votes(1, :)) - 1;
    
    result = zeros(rows, cols);
    
    for i = 1:rows
        currRow = votes(i, :);
        result(i, :) = currRow(currRow ~= losingCandidate);
    end
    
    votes = result;
    
end


%% Discussion Questions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
comment start

1) Console Output

    Vote 1 Data:
                           1     2     3     4     5     6     7     8  
    Round 1 Totals:     1789   518   536   412  1865  1768   287   325  
    Round 2 Totals:     1839   559   579   450  1898  1811     0   364  
    Round 3 Totals:     1901   618   635   507  1958  1881     0     0  
    Round 4 Totals:     2006   713   741     0  2058  1982     0     0  
    Round 5 Totals:     2176     0   924     0  2235  2165     0     0  
    Round 6 Totals:     2472     0     0     0  2534  2494     0     0  
    Round 7 Totals:        0     0     0     0  3737  3763     0     0  
    Winning Candidate:  6


    Vote 2 Data:
                           1     2     3     4     5     6  
    Round 1 Totals:      608   897  1080   936   775   704  
    Round 2 Totals:        0  1015  1212  1047   901   825  
    Round 3 Totals:        0  1225  1426  1245  1104     0  
    Round 4 Totals:        0  1603  1794  1603     0     0  
    Round 5 Totals:        0     0  2600  2400     0     0  
    Winning Candidate:  3

2) Vote 2 Data Modification: Candidate #4 Now Gets Removed During Round 4 Instead of #2

    Vote 2 Data Results (after modification):

                           1     2     3     4     5     6  
    Round 1 Totals:      608   897  1080   936   775   704  
    Round 2 Totals:        0  1015  1212  1047   901   825  
    Round 3 Totals:        0  1225  1426  1245  1104     0  
    Round 4 Totals:        0  1603  1794  1603     0     0  
    Round 5 Totals:        0  2353  2647     0     0     0  
    Winning Candidate:  3

    The difference is not too significant since candidate 3 still wins and
    there are 5 rounds in total. But do note that candidate #2 now survives
    until round 5. 

3) Weighted Ranked-Voting System

    Vote 1 Data:

    Totals:     37124  32025  32550  31442  37389  37337  30691  31442
    Winning Candidate: 5


    Vote 2 Data:

    Totals:     16919  17597  18326  17677  17150  17331
    Winning Candidate: 3

    Note that using the weighted voting system changes the votes1.mat data
    to yield candidate 5 as the winner instead of candidate 6, while the
    winner for votes2.mat stayed the same. The weighted voting system is
    generally better than the normal ranked voting system because the
    ordering of voter preferences have greater importance and influence the
    results of the election to a greater extent than the normal ranked
    version. In the unweighted version, we see exceptional cases where a
    candidate initially in the minority can win the election through many
    rounds of elimination of his/her opponents, which may be unintended by
    the populace. 

comment end
%}
