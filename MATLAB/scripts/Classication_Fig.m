% Region wise
labels_r = [repmat(0,1,33), repmat(1,1,3)]' ;
scores_r = [repmat(0,1,16),repmat(1,1,17),repmat(1,1,3)]'; % Predicted scores

% Subject Wıse
labels_s = [ repmat(1,1,10), repmat(0,1,10) ] ;
scores_s = [1,1,1,1,1,1,1,0,0,0....
    1,1,0,0,0,0,0,0,0,0] ;


% Calculate ROC curve using perfcurve
[falsePositiveRate_s, truePositiveRate_s, ~, AUC_s] = perfcurve(labels_s, scores_s, 1);
[falsePositiveRate_r, truePositiveRate_r, ~, AUC_r] = perfcurve(labels_r, scores_r, 1);

% Plot ROC curve
figure;
plot(falsePositiveRate_s, truePositiveRate_s, 'b-', 'LineWidth', 3);
axis equal;
hold on
plot(falsePositiveRate_r, truePositiveRate_r, 'r-', 'LineWidth', 3);
axis equal;
title('Receiver Operating Characteristic (ROC) Curve','FontSize',20);
xlabel('False Positive Rate','FontSize',30);
ylabel('True Positive Rate','FontSize',30);
grid on;

% Axis option
ax = gca ;
ax.LineWidth = 1 ;
ax.FontSize = 30 ;
ax.XTick = [-1 0 1]; 

% Display AUC value
text(0.4, 0.9, ['AUC = ', num2str(AUC_r)], 'FontSize', 30, 'HorizontalAlignment', 'center','Color','red');
text(0.2, 0.8, ['AUC = ', num2str(AUC_s)], 'FontSize', 30, 'HorizontalAlignment', 'center','Color','blue');

% Optionally display a diagonal line representing random guessing
plot([0, 1], [0, 1], '--', 'Color', [0.7, 0.7, 0.7], 'LineWidth', 1.5);
legend('ROC Curve (Subject-Wise)', 'ROC Curve (Region-Wise)','Random Guessing','Location', 'southeast','FontSize',25);

%text(-0.9, 1, 'C', 'FontSize', 40, 'FontWeight', 'bold');

hold off;

print('E:/EIB/FIGURES/ROC','-dpng','-r300');
