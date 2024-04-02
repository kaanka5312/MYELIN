% Assuming 'dataTable' is your table with the 4th column as the class column ('Class1')

% Calculate initial differences between means for the first three columns
initial_differences = zeros(1,3);
for j = 1:3
    data_array = dataTable{:,j};
    means_1 = mean(data_array(dataTable.Class1 == 1));
    means_0 = mean(data_array(dataTable.Class1 == 2)); % Assuming 0s, based on initial question
    initial_differences(j) = means_1 - means_0;
end

% Bootstrap process
num_bootstrap = 100000;
bootstrap_differences = zeros(num_bootstrap, 3); % Store differences for each bootstrap

% Controlled shuffling within constraints
for t = 1:3
    for i = 1:num_bootstrap
        % Prepare for shuffling the class labels in the 4th column with constraints
        shuffled_classes = dataTable.Class1;
        % Implement a controlled shuffling mechanism here 
        % (Details of this mechanism will follow below)
        % Assuming a function 'controlled_shuffle' that takes the class column and returns a shuffled version
        shuffled_classes = controlled_shuffle(shuffled_classes-1);
        
        % Recalculate means and differences for shuffled data
        data_array = dataTable{:,t};
        means_1_shuffled = mean(data_array(shuffled_classes == 1));
        means_0_shuffled = mean(data_array(shuffled_classes == 0)); % Assuming 0s, adjust if necessary
        bootstrap_differences(i, t) = means_1_shuffled - means_0_shuffled;
    end
end

% Calculate p-values for each column considering the direction
p_values = zeros(1,3);
for j = 1:3
    if initial_differences(j) > 0
        p_values(j) = sum(bootstrap_differences(:,j) <= initial_differences(j)) / num_bootstrap;
    else
        p_values(j) = sum(bootstrap_differences(:,j) >= initial_differences(j)) / num_bootstrap;
    end
end

% Display initial differences and p-values
disp('Initial differences:');
disp(initial_differences);
disp('P-values:');
disp(p_values);
