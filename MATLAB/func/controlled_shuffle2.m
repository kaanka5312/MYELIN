function shuffled_classes = controlled_shuffle2(classes)
    % Extract the indices of 1s
    ones_indices = find(classes == 1);
    
    % Calculate the total number of 1s and the length of the array
    num_ones = length(ones_indices);
    total_length = length(classes);
    
    % Determine the number of 1s for each half
    ones_first_half = floor(num_ones / 2);
    ones_second_half = num_ones - ones_first_half;
    
    % Initialize the shuffled_classes array
    shuffled_classes = zeros(size(classes));
    
    % Initialize counters for how many 1s have been placed in each half
    count_first_half = 0;
    count_second_half = 0;
    
    % Loop to distribute 1s, aiming to balance between first and second half
    for i = 1:num_ones
        % Determine target half based on the current distribution of 1s
        if count_first_half < ones_first_half
            % Target the first half
            block_start = 1;
            block_end = floor(total_length / 2);
            count_first_half = count_first_half + 1;
        else
            % Target the second half
            block_start = floor(total_length / 2) + 1;
            block_end = total_length;
            count_second_half = count_second_half + 1;
        end
        
        % Select a random index within the target block for a 1
        rand_index_within_block = randi([block_start, block_end], 1, 1);
        
        % If the selected position is already taken, find the next available position
        while shuffled_classes(rand_index_within_block) == 1
            rand_index_within_block = rand_index_within_block + 1;
            if rand_index_within_block > block_end % Wrap within the current half
                rand_index_within_block = block_start;
            end
        end
        
        % Place a 1 at the selected position
        shuffled_classes(rand_index_within_block) = 1;
    end
    
    % The rest are implicitly 0s, but you can explicitly set them if needed
    % shuffled_classes(shuffled_classes == 0) = 0; % Optional for clarity
end
