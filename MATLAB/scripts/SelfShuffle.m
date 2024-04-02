% Assume 'matrix' is your 360x4 matrix with the 4th column as the class column

% Calculate initial differences between means for the first three columns

initial_differences = zeros(1,3) ;
for j = 1:3
    data_array = dataTable{:,j} ;
    means_1 = mean(data_array(dataTable.Class1==1)) ;
    means_0 = mean(data_array(dataTable.Class1==2)) ;
    initial_differences(j) = means_1 - means_0 ;
end


% Bootstrap process
num_bootstrap = 10000;
bootstrap_differences = zeros(num_bootstrap, 3); % Store differences for each bootstrap
for t = 1:3
for i = 1:num_bootstrap
    % Shuffle the class labels in the 4th column
    shuffled_classes = dataTable.Class1 ;
    shuffled_classes = shuffled_classes(randperm(length(shuffled_classes)));
    
    % Recalculate means and differences for shuffled data
    data_array = dataTable{:,t} ;
    means_1_shuffled = mean(data_array(shuffled_classes==1)) ;
    means_0_shuffled = mean(data_array(shuffled_classes==2)) ;

    bootstrap_differences(i, t) = means_1_shuffled - means_0_shuffled;
end
end

% Calculate p-values for each column
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
