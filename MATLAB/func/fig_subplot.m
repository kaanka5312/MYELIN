function f = fig_subplot(x, y, col, r, p, xTitle)     

% Fit a linear regression model
lm = fitlm(x, y);

% Plot the scatter plot
f = scatter(x, y, 20, col, 'filled','DisplayName', 'Data');
hold on;

% Plot the regression line
xValues = linspace(min(x), max(x), 100);
yFit = predict(lm, xValues');

if p<0.05
plot(xValues, yFit, 'r-', 'LineWidth', 2, 'DisplayName', strcat('r=', num2str(r) ,' p_{spin} < 0.05'));
else
plot(xValues, yFit, 'r-', 'LineWidth', 2, 'DisplayName', strcat('r=', num2str(r) ,' p_{spin} > 0.05'));
end
% Plot confidence intervals
[yPred, delta] = predict(lm, xValues', 'Prediction', 'curve', 'Confidence', 'on');
plot(xValues, yPred, 'b--', 'DisplayName', '95% Confidence Interval');
plot(xValues, yPred + delta, 'b--', 'HandleVisibility', 'off');
plot(xValues, yPred - delta, 'b--', 'HandleVisibility', 'off');

% Customize the plot
ylabel(xTitle,'FontSize',10);
xlabel('Intracortical Myelin Content','FontSize',10);
%title('Spearman Correlation (p_{spin} < 0.05, r= 0.309)', 'FontSize', 14);
legend('show', 'FontSize', 6);
