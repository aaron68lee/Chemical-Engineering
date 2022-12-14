
function bool = birthdayTest(group)

    % function birthdayTest returns a boolean value
    % true if two people in the vector, group, have birthdays within the same
    % week (any two values are less than 7 from each other), and false
    % otherwise
    
    % same week condition automatically fails if no pair exists
    
    if length(group) <= 1
        bool = false;
        return;
    end
    
    % birthday of last person added
    
    new = group(end);
    
    %fprintf("Last added: %i\n", new);
    
    % loop through all people in the group array and compare their
    % birthdays to the last added person
    
    for i = 1:length(group) - 1
        
        % birthdays must be less than 7 days away
        
        if (abs(group(i) - new) < 7) || (abs(group(i) - (new + 365)) < 7)
            bool = true;
            return;
        end
    end
    
    % no such same-week pair of birthdays was found
    
    bool = false;
    
end