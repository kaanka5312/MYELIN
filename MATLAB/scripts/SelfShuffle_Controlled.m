% Assuming 'dataTable' is your table with the 4th column as the class column ('Class1')
load('C:/Users/kaan/Documents/NatComm2023/MYELIN/DATA/MED_TABLE.mat')

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

% Set up a figure for the plots
p_values = 1 - p_values ;
fig = figure;

% Number of measures to plot
numMeasures = size(bootstrap_differences, 2);

% Plot each measure in a subplot
for i = 1:numMeasures
    subplot(1, numMeasures, i);
    
    % Plot histogram of bootstrap differences
    histogram(bootstrap_differences(:, i), 50, 'FaceColor', 'blue', 'EdgeColor', 'black');
    hold on;
    
    % Plot initial difference as a red dashed line
    xline(initial_differences(i), 'r--', 'LineWidth', 2);
    
    % Title and labels
    %title(sprintf('Measure %d\nInitial Diff: %.2f, p-value: %.4f', i, initial_differences(i), p_values(i)));
     % Get the name of the current measure from the dataTable
    measureName = dataTable.Properties.VariableNames{i};
    
    % Title and labels using the measure name
    title(sprintf('%s\nInitial Diff: %.2f, p-value: %.4f', measureName, initial_differences(i), p_values(i)));
    xlabel('Bootstrap Difference');
    if i == 1
        ylabel('Frequency');
    end
 
end

% Adjust subplot layouts to avoid overlap
sgtitle('Controlled Bootstrapping and Initial Empirical Values');

% Save the figure with more options
exportgraphics(fig, 'C:/Users/kaan/Desktop/Controlled_Shuffle.png', 'Resolution', 300);

