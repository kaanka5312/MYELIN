% Reading each tseries file
mainFolder = '/media/kaansocat/Elements/EIB/';

% Get the directory listing
contents = dir(mainFolder);

% Initialize an empty cell array to store folder names
folderNames = {};

% Iterate through each item in the directory
for i = 1:numel(contents)
    % Check if the item is a directory (folder)
    if contents(i).isdir && ~strcmp(contents(i).name, '.') && ~strcmp(contents(i).name, '..')
        % Add the folder name to the cell array
        folderNames = [folderNames; contents(i).name];
    end
end

% Each subject
for i=1:100
CIFTI = cifti_read([mainFolder,folderNames{i},'/INT/Glasser_cortical_timeseries.ptseries.nii']) ;

ACW50_arr = [] ; 
ACW0_arr = [] ;
% Calculating ACW-50 and ACW-0 for each parcellation 
for t=1:360
    [ACW0, ACW50] = ACW_kaan(CIFTI.cdata(t,:),1/0.72,100,50) ;

    ACW50_arr(t) = ACW50 ;
    ACW0_arr(t) = ACW0 ;
end
% Saving each subject as a row parcellation as columns 

save([mainFolder,folderNames{i},'/INT/INT.mat'], 'ACW0_arr' , 'ACW50_arr');
end

%%%====== AVERAGING FOR SPINNING =========

% ACW data 
ACW50_all = [] ;
ACW0_all = [] ;

for i=1:100

load([mainFolder,folderNames{i},'/INT/INT.mat']) ;
ACW50_all = vertcat(ACW50_all, ACW50_arr) ;
ACW0_all = vertcat(ACW0_all, ACW0_arr) ;

end

% Myelin data 

myelin_all = [] ;

for i=1:100

MY = load([mainFolder,folderNames{i},'/Glasser_myelin_content.txt']) ;
myelin_all = vertcat(myelin_all, MY') ;

end

save([mainFolder,'INT_all.mat'], 'ACW0_all', 'ACW50_all', 'myelin_all')


% Spin permutation testing for two cortical maps
AV_ACW50 = mean(ACW50_all,1) ;
AV_ACW0 = mean(ACW0_all,1) ;
AV_MY = mean(myelin_all,1) ;


% Git URL for the function for spin test. LATER VERIFY WITH ENIGMA TOOLBOX
% https://github.com/frantisekvasa/rotate_parcellation.git
% Coordinates for glasser parcellations (HCP)
CORD = load('/home/kaansocat/Documents/MATLAB/rotate_parcellation-master/sphere_HCP.txt') ;
perm_id = rotate_parcellation(CORD(1:180,:), CORD(181:360,:), 10.000) ;

% Generates p valie from above 
perm_sphere_p(AV_ACW0', AV_MY', perm_id, 'spearman')

% ACW0 content to myelin content
% Example data
x = AV_ACW0' ;
y = AV_MY' ;

% See the r-value
corr(x,y,'type','Spearman')

% Fit a linear regression model
lm = fitlm(x, y);

% Plot the scatter plot
scatter(x, y, 'filled', 'DisplayName', 'Data');
hold on;

% Plot the regression line
xValues = linspace(min(x), max(x), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'r-', 'LineWidth', 2, 'DisplayName', 'Regression Line');

% Plot confidence intervals
[yPred, delta] = predict(lm, xValues', 'Prediction', 'curve', 'Confidence', 'on');
plot(xValues, yPred, 'b--', 'DisplayName', '95% Confidence Interval');
plot(xValues, yPred + delta, 'b--', 'HandleVisibility', 'off');
plot(xValues, yPred - delta, 'b--', 'HandleVisibility', 'off');

% Customize the plot
xlabel('ACW-0 (seconds)','FontSize',20);
ylabel('Intracortical Myelin Content','FontSize',20);
title('Spearman Correlation (p_{spin} < 0.001, r= -0.68)', 'FontSize', 24);
legend('show');
grid on;
hold off;


% ACW50 content to myelin content
% Example data
x = AV_ACW50' ;
y = AV_MY' ;

% Fit a linear regression model
lm = fitlm(x, y);

% See the r-value
corr(x,y,'type','Spearman')

% Plot the scatter plot
scatter(x, y, 'filled', 'DisplayName', 'Data');
hold on;

% Plot the regression line
xValues = linspace(min(x), max(x), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'r-', 'LineWidth', 2, 'DisplayName', 'Regression Line');

% Plot confidence intervals
[yPred, delta] = predict(lm, xValues', 'Prediction', 'curve', 'Confidence', 'on');
plot(xValues, yPred, 'b--', 'DisplayName', '95% Confidence Interval');
plot(xValues, yPred + delta, 'b--', 'HandleVisibility', 'off');
plot(xValues, yPred - delta, 'b--', 'HandleVisibility', 'off');

% Customize the plot
xlabel('ACW-0 (seconds)','FontSize',20);
ylabel('Intracortical Myelin Content','FontSize',20);
title('Spearman Correlation (p_{spin} = 0.20, r= 0.1681)', 'FontSize', 24);
legend('show');
grid on;
hold off;

% ============== SELF REGIONS ===============
T = readtable('C:/Users/kaan/Downloads/Self2Glasser.csv') ;

INT = sort([108, 220, 120, 302, 286, 291, 148, 63, 258, 189]) ;
EXT = sort([138, 82, 78, 109, 318, 145, 116, 48, 292, 20, 297, 230, 57, 249]) ;
MENT = sort([241, 214, 291, 330, 248, 212, 312, 150, 107, 267, 78]) ;


%%
% Self Regions that each group is inbetween
CORD = load('C:/Users/kaan/Documents/MATLAB/rotate_parcellation-master/rotate_parcellation-master/sphere_HCP.txt') ;

SELF = sort([INT,EXT,MENT]) ;

% Droppin one of double 78 
indexToRemove = 5;

% Remove the element at the specified index
SELF = SELF([1:indexToRemove-1, indexToRemove+1:end]);

% Self regions as a whole, without dividing according to layers 

perm_id = rotate_parcellation(CORD(SELF(1:16),:), CORD(SELF(17:34),:), 10.000) ;
perm_id = rotate_parcellation(CORD(INT(1:4),:), CORD(INT(5:10),:), 100.000) ;

% Generates p value from above 
perm_sphere_p(mean(ACW0_all(:,SELF),1)', mean(myelin_all(:,SELF),1)', perm_id, 'spearman')

% ACW0 content to myelin content
% Example data
x = AV_ACW0' ;
y = AV_MY' ;

% See the r-value
corr(x,y,'type','Spearman')

% Fit a linear regression model
lm = fitlm(x, y);

% Plot the scatter plot
scatter(x, y, 'filled', 'DisplayName', 'Data');
hold on;

% Plot the regression line
xValues = linspace(min(x), max(x), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'r-', 'LineWidth', 2, 'DisplayName', 'Regression Line');

% Plot confidence intervals
[yPred, delta] = predict(lm, xValues', 'Prediction', 'curve', 'Confidence', 'on');
plot(xValues, yPred, 'b--', 'DisplayName', '95% Confidence Interval');
plot(xValues, yPred + delta, 'b--', 'HandleVisibility', 'off');
plot(xValues, yPred - delta, 'b--', 'HandleVisibility', 'off');

% Customize the plot
xlabel('ACW-0 (seconds)','FontSize',20);
ylabel('Intracortical Myelin Content','FontSize',20);
title('Spearman Correlation (p_{spin} < 0.001, r= -0.68)', 'FontSize', 24);
legend('show');
grid on;
hold off;




% Interoceptive 
perm_id = rotate_parcellation(CORD(INT(1:4),:), CORD(INT(5:10),:), 100.000) ;
% Generates p valie from above 
perm_sphere_p(mean(ACW0_all(:,INT),1)', mean(myelin_all(:,INT),1)', perm_id, 'spearman')

% Exteroceptive 
perm_id = rotate_parcellation(CORD(EXT(1:9),:), CORD(EXT(10:14),:), 100.000) ;
% Generates p valie from above 
perm_sphere_p(mean(ACW0_all(:,EXT),1)', mean(myelin_all(:,EXT),1)', perm_id, 'spearman')

% Mental 
perm_id = rotate_parcellation(CORD(MENT(1:3),:), CORD(MENT(4:11),:), 100.000) ;
% Generates p valie from above 
perm_sphere_p(mean(ACW0_all(:,MENT),1)', mean(myelin_all(:,MENT),1)', perm_id, 'spearman')

%%============= FIGURE OF SELF REGIONS ACW-0 to Myelin Content ============
scatter(mean(ACW0_all(:,INT),1), mean(myelin_all(:,INT),1), 'filled','red','DisplayName','Interoceptive' )
hold on 

% Fit a linear regression model
lm = fitlm(mean(ACW0_all(:,INT),1), mean(myelin_all(:,INT),1));
% Plot the regression line
xValues = linspace(min(mean(ACW0_all(:,INT),1)), max(mean(ACW0_all(:,INT),1)), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'r-', 'LineWidth', 2, 'DisplayName', 'r_{INT}= -0.70, p_{spin} < 0.05');
corr(mean(ACW0_all(:,INT),1)', mean(myelin_all(:,INT),1)','type','Spearman' );

scatter(mean(ACW0_all(:,EXT),1), mean(myelin_all(:,EXT),1), 'filled','blue','DisplayName','Exteroceptive'  )

% Fit a linear regression model
lm = fitlm(mean(ACW0_all(:,EXT),1), mean(myelin_all(:,EXT),1));
% Plot the regression line
xValues = linspace(min(mean(ACW0_all(:,EXT),1)), max(mean(ACW0_all(:,EXT),1)), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'b-', 'LineWidth', 2, 'DisplayName', 'r_{EXT}= -0.48, p_{spin} > 0.05');
[RHO,PVAL]=corr(mean(ACW0_all(:,EXT),1)', mean(myelin_all(:,EXT),1)','type','Spearman' ) ;

scatter(mean(ACW0_all(:,MENT),1), mean(myelin_all(:,MENT),1), 'filled','greed', 'DisplayName','Mental' )

% Fit a linear regression model
lm = fitlm(mean(ACW0_all(:,MENT),1), mean(myelin_all(:,MENT),1));
% Plot the regression line
xValues = linspace(min(mean(ACW0_all(:,MENT),1)), max(mean(ACW0_all(:,MENT),1)), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'g-', 'LineWidth', 2, 'DisplayName', 'r_{MENT}= -0.77, p_{spin} < 0.05');
[RHO,PVAL]=corr(mean(ACW0_all(:,MENT),1)', mean(myelin_all(:,MENT),1)','type','Spearman' );

% Customize the plot
xlabel('ACW-0 (seconds)','FontSize',20);
ylabel('Intracortical Myelin Content','FontSize',20);
title('Spearman Correlation', 'FontSize', 24);
legend('show','FontSize', 20);
grid on;
hold off

%===================== END OF SELF THREE LAYER CORR FIGURE ================== 





% Example data with three groups
group1 = reshape(myelin_all(:,INT),1,[])';
group2 = reshape(myelin_all(:,EXT),1,[])';
group3 = reshape(myelin_all(:,MENT),1,[])';

% Store data in a cell array
data = {group1, group2, group3};

% Find the maximum length among the vectors in the cell array
maxLen = max(cellfun(@length, data));

% Pad vectors with NaN to make them of equal length
dataPadded = cellfun(@(x) [x; nan(maxLen - length(x), 1)], data, 'UniformOutput', false);

% Convert the padded cell array to a matrix
dataMatrix = cell2mat(dataPadded);

% Create a boxplot
figure;
boxplot(dataMatrix, 'Labels', {'Group 1', 'Group 2', 'Group 3'});
title('Boxplot of Three Groups');

% Perform one-way ANOVA with unequal sample sizes
[pValue, ~, stats] = anova1(dataMatrix, [], 'on');

% Display ANOVA results
fprintf('One-Way ANOVA p-value: %f\n', pValue);

% If ANOVA is significant, perform post hoc comparisons using Bonferroni method
if pValue < 0.05
    comparisonResults = multcompare(stats, 'CType', 'bonferroni');
    fprintf('Post Hoc Comparisons:\n');
    disp(comparisonResults);
    
    % Compare means
    groupMeans = grpstats(dataMatrix, ones(size(dataMatrix)), 'mean');
    fprintf('Group Means:\n');
    disp(groupMeans);
else
    disp('No significant difference detected.');
end
