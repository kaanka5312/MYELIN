function shuffled_classes = controlled_shuffle(classes)
    % Extract the indices of 1s and 0s
    ones_indices = find(classes == 1);
    zeros_indices = find(classes == 0);
    
    % Define blocks to ensure no clustering at the extremes
    num_ones = length(ones_indices);
    total_length = length(classes);
    block_size = floor(total_length / num_ones); % Adjust block size based on your needs
    
    % Initialize the shuffled_classes array
    shuffled_classes = zeros(size(classes));
    
    % Distribute 1s across the array to avoid clustering
    for i = 1:num_ones
        % Calculate block start and end indices
        block_start = (i-1)*block_size + 1;
        block_end = i*block_size;
        
        % Ensure we do not exceed the array bounds
        if block_start > total_length
            block_start = total_length;
        end
        if block_end > total_length
            block_end = total_length;
        end
        
        % Select a random index within the block for a 1
        rand_index_within_block = randi([block_start, block_end], 1, 1);
        
        % If the selected position is already taken, find the next available position
        while shuffled_classes(rand_index_within_block) == 1
            rand_index_within_block = rand_index_within_block + 1;
            if rand_index_within_block > total_length
                rand_index_within_block = 1; % Restart from the beginning if we reach the end
            end
        end
        
        % Place a 1 at the selected position
        shuffled_classes(rand_index_within_block) = 1;
    end
    
    % Fill the rest with 0s (if any index remains unset, this step ensures they are set to 0)
    shuffled_classes(shuffled_classes == 0) = 0; % This line may seem redundant but clarifies intent
    
end




