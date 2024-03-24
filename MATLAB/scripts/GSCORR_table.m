%% This script is only important for it produces MED file in the end
% save('E:/EIB/MED.mat','MED'). MED can be find in DATA folder
% Author: kaanka5312
%-----%-----%-----%-----%-----%-----%-----%-----%-----%-----

% Reading each tseries file
mainFolder = 'E:/EIB/100_Subj/';

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


GSCORR_arr = zeros(360,100) ;

% Each subject
for i=1:100
    % Importing whole graymatter data 
    GRAYMATTER = cifti_read([mainFolder,folderNames{i},'/INT/avg_norm.dtseries.nii']) ;
    GS = mean(GRAYMATTER.cdata,1) ; 
    GLASSER = cifti_read([mainFolder,folderNames{i},'/INT/Glasser_cortical_timeseries.ptseries.nii']) ;
    GSCORR_arr(:,i) = GSCORR(GS, GLASSER.cdata) ;
end
% Saving each subject as a row parcellation as columns 

save([mainFolder,folderNames{i},'/INT/GSCORR.mat'], 'GSCORR_arr');

%% Whole Brain
load('E:/EIB/100_Subj/GSCORR.mat')
load('E:/EIB/MATLAB_ClassificationApp/INT_all.mat')

AV_MY = mean(myelin_all,1) ;
AV_GSCORR = mean(GSCORR_arr,2) ;


CORD = load('C:/Users/kaan/Documents/MATLAB/rotate_parcellation-master/rotate_parcellation-master/sphere_HCP.txt') ;
perm_id = rotate_parcellation(CORD(1:180,:), CORD(181:360,:), 10000 ) ;

% Generates p valie from above 
[p_spin, r_dist] = perm_sphere_p(AV_GSCORR, AV_MY', perm_id, 'spearman') ;

% Store p-values and null distributions
p_and_d =  cell2struct({[p_spin; r_dist]}, ...
                       {'GSCORR_MY'}, 2);

f = figure,
    set(gcf,'color','w');
    set(gcf,'units','normalized','position',[0 0 1 0.3])
    fns = fieldnames(p_and_d);

    for k = 1:numel(fieldnames(rvals))
        % Define plot colors
        if k <= 2; col = [0.66 0.13 0.11]; else; col = [0.2 0.33 0.49]; end

        % Plot null distributions
        axs = subplot(1, 4, k); hold on
        h = histogram(p_and_d.(fns{k})(2:end), 50, 'Normalization', 'pdf', 'edgecolor', 'w', ...
                      'facecolor', col, 'facealpha', 1, 'linewidth', 0.5);
        l = line([rvals.(fns{k}) rvals.(fns{k})], get(gca, 'ylim'), 'linestyle', '--', ...
                 'color', 'k', 'linewidth', 1.5);
        xlabel(['Null correlations' newline '(' strrep(fns{k}, '_', ' ') ')'])
        ylabel('Density')
        legend(l,['{\it r}=' num2str(round(rvals.(fns{k}), 2)) newline ...
                  '{\it p}=' num2str(round(p_and_d.(fns{k})(1), 3))])
        legend boxoff
    end



% GSCORR to myelin content
% Example data
x = AV_GSCORR ;
y = AV_MY' ;

% See the r-value
[RHO,PVAL]=corr(x,y,'type','Spearman')

% Fit a linear regression model
lm = fitlm(x, y);

% Plot the scatter plot
scatter(x, y, 'filled', 'black', 'DisplayName', 'Data');
hold on;

% Plot the regression line
xValues = linspace(min(x), max(x), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'r-', 'LineWidth', 2, 'DisplayName', 'r=0.309, p_{spin} < 0.05');

% Plot confidence intervals
[yPred, delta] = predict(lm, xValues', 'Prediction', 'curve', 'Confidence', 'on');
plot(xValues, yPred, 'b--', 'DisplayName', '95% Confidence Interval');
plot(xValues, yPred + delta, 'b--', 'HandleVisibility', 'off');
plot(xValues, yPred - delta, 'b--', 'HandleVisibility', 'off');

% Customize the plot
xlabel('GSCORR','FontSize',20);
ylabel('Intracortical Myelin Content','FontSize',20);
title('Spearman Correlation (p_{spin} < 0.05, r= 0.309)', 'FontSize', 24);
legend('show', 'FontSize', 20);
grid on;
hold off;

savefig('E:/EIB/FIGURES/Global.fig')
%% All self regions with GSCORR 

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

% Droppin one of double 291 
indexToRemove = 27;

% Remove the element at the specified index
SELF = SELF([1:indexToRemove-1, indexToRemove+1:end]);
% Self regions as a whole, without dividing according to layers 

% Spin permutation testing for two cortical maps
perm_id = rotate_parcellation(CORD(SELF(1:15),:), CORD(SELF(16:end),:), 100.000) ;
% perm_id = rotate_parcellation(CORD(INT(1:4),:), CORD(INT(5:10),:), 100.000) ;

% Generates p value from above 
perm_sphere_p(AV_GSCORR(SELF), AV_MY(SELF)', perm_id, 'spearman')

% GSCORR to myelin content
x = AV_GSCORR(SELF) ;
y = AV_MY(SELF)' ;

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
xlabel('GSCORR','FontSize',20);
ylabel('Intracortical Myelin Content','FontSize',20);
title('Three Layer of Self (p_{spin} < 0.05, r=0.52)', 'FontSize', 24);
legend('show');
grid on;
hold off;

savefig('E:/EIB/FIGURES/Self.fig')

% Non Self Regions
NONSELF = 1:360 ;
NONSELF = setdiff( NONSELF, SELF) ;

perm_id = rotate_parcellation(CORD(NONSELF(1:165),:), CORD(NONSELF(166:end),:), 100.000) ;

% Generates p value from above 
perm_sphere_p(AV_GSCORR(NONSELF), AV_MY(NONSELF)', perm_id, 'spearman')

x = AV_GSCORR(NONSELF) ;
y = AV_MY(NONSELF)' ;

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
xlabel('GSCORR','FontSize',20);
ylabel('Intracortical Myelin Content','FontSize',20);
title('Non-Self (p_{spin} > 0.05, r=0.29)', 'FontSize', 24);
legend('show');
grid on;
hold off;

savefig('E:/EIB/FIGURES/Nonself.fig')


h1 = openfig('E:/EIB/FIGURES/Global.fig','reuse'); % open figure
ax1 = gca; % get handle to axes of figure
h2 = openfig('E:/EIB/FIGURES/Self.fig','reuse'); % open figure
ax2 = gca; % get handle to axes of figure
h3 = openfig('E:/EIB/FIGURES/Nonself.fig','reuse');
ax3 = gca;
% test1.fig and test2.fig are the names of the figure files which you would % like to copy into multiple subplots

h4 = figure; %create new figure
s1 = subplot(1,3,1); %create and get handle to the subplot axes
title('Global (p_{spin} < 0.05, r= 0.309)', 'FontSize', 16);
ylabel('Intracortical Myelin Content','FontSize', 20)
s2 = subplot(1,3,2); %create and get handle to the subplot axes
title('Three Layer of Self (p_{spin} < 0.05, r=0.52)', 'FontSize', 16);
xlabel('GSCORR', 'FontSize', 20)
s3 = subplot(1,3,3);
title('Non-Self (p_{spin} > 0.05, r=0.29)', 'FontSize', 16);
fig1 = get(ax1,'children'); %get handle to all the children in the figure
fig2 = get(ax2,'children');
fig3 = get(ax3,'children');

copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
copyobj(fig2,s2);
copyobj(fig3,s3); %copy children to new parent axes i.e. the subplot axes

%% Each Layer Respectively

% Alternatively, include intersection regions to only one layer
% Including 291 to only interoceptive
INT = sort([108, 220, 120, 302, 286, 291, 148, 63, 258, 189]) ;
EXT = sort([138, 82, 78, 109, 318, 145, 116, 48, 292, 20, 297, 230, 57, 249]) ;
% Including 78 to only exteroceptive 
MENT = sort([241, 214, 330, 248, 212, 312, 150, 107, 267]) ;

% Interoceptive 
perm_id = rotate_parcellation(CORD(INT(1:4),:), CORD(INT(5:10),:), 100.000) ;
% Generates p valie from above 
perm_sphere_p(AV_GSCORR(INT), AV_MY(INT)', perm_id, 'spearman')

% Exteroceptive 
perm_id = rotate_parcellation(CORD(EXT(1:9),:), CORD(EXT(10:14),:), 100.000) ;
% Generates p valie from above 
perm_sphere_p(AV_GSCORR(EXT), AV_MY(EXT)', perm_id, 'spearman')

% Mental 
% perm_id = rotate_parcellation(CORD(MENT(1:3),:), CORD(MENT(4:11),:), 100.000) ;
% Generates p valie from above 
% perm_sphere_p(AV_GSCORR(MENT), AV_MY(MENT)', perm_id, 'spearman')

% Mental Alternative. Intersected regions included to lower levels. 
perm_id = rotate_parcellation(CORD(MENT(1:2),:), CORD(MENT(3:end),:), 100.000) ;
% Generates p valie from above 
perm_sphere_p(AV_GSCORR(MENT), AV_MY(MENT)', perm_id, 'spearman')

%%============= FIGURE OF SELF REGIONS ACW-0 to Myelin Content ============
scatter(AV_GSCORR(INT), AV_MY(INT)', 'filled','red','DisplayName','Interoceptive' )
hold on 

% Fit a linear regression model
lm = fitlm( AV_GSCORR(INT), AV_MY(INT)' );
% Plot the regression line
xValues = linspace(min(AV_GSCORR(INT)), max(AV_GSCORR(INT)), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'r-', 'LineWidth', 2, 'DisplayName', 'r_{INT}= 0.41, p_{spin} > 0.05');
corr(AV_GSCORR(INT), AV_MY(INT)','type','Spearman' );

scatter(AV_GSCORR(EXT), AV_MY(EXT)', 'filled','blue','DisplayName','Interoceptive' )
hold on 

% Fit a linear regression model
lm = fitlm( AV_GSCORR(EXT), AV_MY(EXT)' );
% Plot the regression line
xValues = linspace(min(AV_GSCORR(EXT)), max(AV_GSCORR(EXT)), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'b-', 'LineWidth', 2, 'DisplayName', 'r_{EXT}= 0.60, p_{spin} < 0.05');
corr(AV_GSCORR(EXT), AV_MY(EXT)','type','Spearman' );

scatter(AV_GSCORR(MENT), AV_MY(MENT)', 'filled','greed', 'DisplayName','Mental' )

% Fit a linear regression model
lm = fitlm( AV_GSCORR(MENT), AV_MY(MENT)' );
% Plot the regression line
xValues = linspace(min(AV_GSCORR(MENT)), max(AV_GSCORR(MENT)), 100);
yFit = predict(lm, xValues');
plot(xValues, yFit, 'g-', 'LineWidth', 2, 'DisplayName', 'r_{MENT}= 0.41, p_{spin} > 0.05');
[RHO,PVAL]=corr(AV_GSCORR(MENT), AV_MY(MENT)' ,'type','Spearman' )

% Customize the plot
xlabel('GSCORR','FontSize',20);
ylabel('Intracortical Myelin Content','FontSize',20);
title('Each Layer Respectively', 'FontSize', 24);
legend('show','FontSize', 20);
grid on;
hold off
%% GSCORR comparison
% Example data with three groups
group1 = reshape(GSCORR_arr(INT,:),1,[])';
group2 = reshape(GSCORR_arr(EXT,:),1,[])';
group3 = reshape(GSCORR_arr(MENT,:),1,[])';

% group1 = AV_GSCORR(INT) ;
% group2 = AV_GSCORR(EXT) ;
% group3 = AV_GSCORR(MENT);

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
boxplot(dataMatrix, 'Labels', {'INT', 'EXT', 'MENT'});
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

    % Change y-axis tick labels
    newLabels = {'MENTAL','EXT','INT'};
    set(gca, 'YTickLabel', newLabels);

    % Add title and x-axis label to the comparison figure
    compFig = gcf;
    compFig.Name = 'Multiple Comparison Figure';
    title('Post Hoc Comparisons (Bonferroni corrected)');
    xlabel('GSCORR');
    ylabel('Three Topography of Self')
else
    disp('No significant difference detected.');
end

%% combining data for mediation investigation 

load('E:/EIB/INT_all.mat')
% Spin permutation testing for two cortical maps
AV_ACW50 = mean(ACW50_all,1) ;
AV_ACW0 = mean(ACW0_all,1) ;
AV_MY = mean(myelin_all,1) ;


MED = horzcat(AV_GSCORR, AV_ACW0', AV_MY', AV_HURST) ;
columnNames = {'GSCORR', 'ACW0', 'MYELIN', 'HURST'};
% Convert matrix to table with column names
% dataTable = array2table(MED, 'VariableNames', columnNames);

outputArray = repmat(1, 1, 360); % Non self 
outputArray(SELF) = 2; % Self

outputArray_2 = repmat(1, 1, 360); % Non self 
outputArray_2(INT) = 2; % Interoceptive
outputArray_2(EXT) = 3; % Interoceptive
outputArray_2(MENT) = 4; % Interoceptive

MED = horzcat(MED, outputArray', outputArray_2'); 
columnNames = {'GSCORR', 'ACW', 'Myelin','Hurst', 'Class1', 'Class2'} ;
dataTable = array2table(MED, 'VariableNames', columnNames);

save('E:/EIB/MED.mat','MED')
save('E:/EIB/MATLAB_ClassificationApp/MED_TABLE.mat','MED')