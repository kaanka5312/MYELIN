function plotEffectSize(group1, group2, TITLE)
    % Inputs
    % group1: Data for group 1
    % group2: Data for group 2

    % Define labels
    labels = {'Self', 'Non-Self'};
    data_labels = [repmat(labels(1), length(group1), 1); repmat(labels(2), length(group2), 1)];

    % Create table for plotting
    plot_dat = table([group1; group2], data_labels, 'VariableNames', {'values', 'class'});

    % Create categorical data for plotting
    x = categorical(plot_dat.class, labels, 'Ordinal', true);
    y = plot_dat.values;
    
    % Define colors
    c = [repmat(hex2rgb("E64B35"), length(group1), 1); repmat(hex2rgb("4DBBD5"), length(group2), 1)];

    % Plot swarmchart and boxchart
    %figure;
    swarmchart(x, y, 100, c, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
    hold on;
    boxchart(x, plot_dat.values);

    % Perform t-test
    [h, p, ~, stats] = ttest(group1, group2);

    % Annotate significance
    %{
    if h
        text(1.5, max([mean(group1), mean(group2)]) + 0.2, '*', 'FontSize', 30, 'HorizontalAlignment', 'center');
    else
        text(1.5, max([mean(group1), mean(group2)]) + 0.2, 'ns', 'FontSize', 20, 'HorizontalAlignment', 'center');
    end
    %}

    if h
        text(1.5, 0.2 , '*', 'FontSize', 30, 'HorizontalAlignment', 'center');
    else
        text(1.5, 0.2 , 'ns', 'FontSize', 20, 'HorizontalAlignment', 'center');
    end

    % Calculate Effect Size
    mean1 = mean(group1);
    mean2 = mean(group2);
    std1 = std(group1);
    std2 = std(group2);
    n1 = length(group1);
    n2 = length(group2);
    pooledStd = sqrt(((n1 - 1) * std1^2 + (n2 - 1) * std2^2) / (n1 + n2 - 2));
    effectSize = (mean1 - mean2) / pooledStd;

    disp(['Cohen''s d: ', num2str(effectSize)]);

    % Calculate y-axis limits with padding
    %{
    yMin = min(y);
    yMax = max(y);
    yRange = yMax - yMin;
    padding = 0.1 * yRange; % 10% padding
    ylim([yMin - padding, yMax + padding]);
    %}
    ylim([-0.1, 0.3])
    modifiedString = strrep(TITLE, '_', ' ');
    
    % Set plot properties
    set(gca, 'FontSize', 15);  % Adjust the font size as needed
    ylabel('GSCORR', 'FontSize', 15);
    
    title(modifiedString, 'FontSize', 15);
    subtitle(['Cohen''s d: ', num2str(effectSize)], 'FontSize', 15);
    hold off;
end

