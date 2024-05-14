% This script produces figure 5 in myelin project. It uses functions for
% Wilcox-Cowan and Balloon-Windkessel simulations.
% Author: kaanka5312

%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----
%% Making self regions myelin less than nonself regions
addpath('C:/Users/kaan/Documents/NatComm2023/MYELIN/MATLAB/func/') ;
W = load('C:/Users/kaan/Documents/NatComm2023/MYELIN/DATA/averageConnectivity_Fpt.mat'); %DTI common matrix
parcelID = W.parcelIDs ;
% Being sure that matrix parcel organization as same
myDataParcels = cell(360,3) ;
%{
for i=1:360
myDataParcels{i,1} = GLASSER.diminfo{1,1}.parcels(i).name ; %GLASSER is from my data structure.
end
% Vectorized version
% myDataParcels{:,1} = arrayfun(@(i) GLASSER.diminfo{1,1}.parcels(i).name, 1:360, 'UniformOutput', false);
dat=load('E:/EIB/MATLAB_ClassificationApp/MED.mat') ; 
dat.MED(:,5) %Self / Non-Self Data

for i=1:360
myDataParcels{i,2} = dat.MED(i,5) ; %GLASSER is from my data structure.
end % 2 is self 1 is nonself

for i=1:360
myDataParcels{i,3} = dat.MED(i,3) ; % Myelin data
end % 2 is self 1 is nonself

%save("E:/EIB/MATLAB_ClassificationApp/myDataParcels.mat",'myDataParcels') % For later
%}
load("C:/Users/kaan/Documents/NatComm2023/MYELIN/DATA/myDataParcels.mat",'myDataParcels') % For later

% Reordering self / non self parcels for the organization of DTI Matrix

DtiMyelin = zeros(360,2) ; %For DTI common matrix Self / Nonself indices
for i=1:360
index = find(strcmp(myDataParcels, [parcelID{i},'_ROI'])) ;
DtiMyelin(i,1) = myDataParcels{index, 2} ;
DtiMyelin(i,2) = myDataParcels{index, 3} ;
end

%% Parameters
b = -2; % Tweat this to have x (firing rate) around 0.4
tau = 0.1;
C = 1;
s = 0; % rest
k = 1;
tspan = 300;
dt = 10 / 1000;


W = 10.^(W.Fpt); % Matrix Log10 olceginde. Kullanmadan once 10.^ yapmak gerekiyor
n = length(W);
W(logical(eye(n))) = 1; % Makes diagonal elemenst equal to 1. 
W(~logical(eye(n))) = (2*360)*(W(~logical(eye(n))) ./ sum(W(~logical(eye(n))), 'all'));

nsims = 100;

sim_fmri_subj = zeros(360, 18002, nsims);
sim_rate_subj = zeros(360, 20002, nsims);

GSCORR_x_normal= zeros(360,nsims);
GSCORR_fmri_normal = zeros(360,nsims);
GSCORR_fmri_hrf_normal = zeros(360,nsims);

ACW_x_normal = zeros(360,nsims);
ACW_bw_normal = zeros(360,nsims);
ACW_hrf_normal = zeros(360,nsims);
  
parfor i = 1:nsims
        tic
        %% Set up W
        %% run simulations
        [x,time] = wilsoncowan_RK2(tau, b, W, k, s, C, tspan);
        sim_fmri = bw(x,200,1/100) ;
       % sim_fmri_hrf = hrfconv(x, time, 10/1000)
        sim_rate_subj(:,:,i) = x ;
        sim_fmri_subj(:,:,i) = sim_fmri ;
        toc
        fprintf("\ni = %d\n", i)
end

% Define sampling frequency and frequencies for bandpass
TR = 0.72; % TR of HCP data
Fs = 1/TR;  % Sampling frequency in Hz, TR is your repetition time in seconds
lowFreq = 0.01;  % Low cutoff frequency in Hz
highFreq = 1 / (2 * TR);  % High cutoff frequency in Hz.  Up to Nyquist Frequency

f1 = 0.01 ;
f2 = 0.69 ;
Wp1 = f1 / (Fs / 2) ; % Normalizing with Nyquist before filter application
Wp2 = f2 / (Fs / 2) ; % Normalizing with Nyquist before filter application
Wp = [Wp1 Wp2];  % Passband: 
Rp = 0.1 ;     % 0.1 dB passband ripple
N = 4;      % Filter order

[b, a] = cheby1(N, Rp, Wp, 'bandpass');

% ACW0
parfor i=1:nsims
    for t=1:360
       % ACW with firing rate
       [ACW0, ~, ~, ~] = acw(sim_rate_subj(t,:,i),100,false);
       ACW_x_normal(t,i) = ACW0 ;
    end
end

parfor i=1:nsims
    % Downsampling fmri simulation to match with HCP data. TR = 2 sampling
    % is sampling the 0.5 Hz. Thus our TR = 0.72, thus we will sampling at
    % Fs (1/Tr). dt = 100 in simulation, so .72 x 100 = 72 should be
    % downsampling rate.
    sim_fmri_downsampled = downsample(sim_fmri_subj(:,:,i)', 72)' ;
    % Detrending the data 
    sim_fmri_downsampled_dt = detrend(sim_fmri_downsampled) ;
    % Bandpassing the data
    filteredData = zeros(size(sim_fmri_downsampled_dt));
    % Apply the filter to each ROI (each voxel/ROI)
    for t = 1:size(sim_fmri_downsampled_dt, 1)
        filteredData(t, :) = filtfilt(b, a, sim_fmri_downsampled_dt(t, :));
    end
    
    % Calculating GS from filtered Data
    GS_fmri = mean(sim_fmri_downsampled_dt,1) ;
    GS_fmri_filtered = filtfilt(b, a, GS_fmri);
    res_corr = corr(filteredData',GS_fmri_filtered') ;
    %Fisher Z transformaiton
    GSCORR_fmri_normal(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); 
end

parfor i=1:nsims
    % Downsampling fmri simulation to match with HCP data. TR = 2 sampling
    % is sampling the 0.5 Hz. Thus our TR = 0.72, thus we will sampling at
    % Fs (1/Tr). dt = 100 in simulation, so .72 x 100 = 72 should be
    % downsampling rate.
    sim_fmri_downsampled = downsample(sim_fmri_subj(:,:,i)', 72)' ;
    % Detrending the data 
    %sim_fmri_downsampled_dt = detrend(sim_fmri_downsampled) ;
    % Bandpassing the data
    %filteredData = zeros(size(sim_fmri_downsampled_dt));
    % Apply the filter to each ROI (each voxel/ROI)
    %for t = 1:size(sim_fmri_downsampled_dt, 1)
     %   filteredData(t, :) = filtfilt(b, a, sim_fmri_downsampled_dt(t, :));
    %end
    
    % Calculating GS from filtered Data
    %GS_fmri = mean(sim_fmri_downsampled_dt,1) ;
    %GS_fmri_filtered = filtfilt(b, a, GS_fmri);
    %res_corr = corr(filteredData',GS_fmri_filtered') ;
    %Fisher Z transformaiton
    %GSCORR_fmri_normal(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); 

    % Calculating GS from filtered Data
    GS_fmri = mean(sim_fmri_downsampled,1) ;
    res_corr = corr(sim_fmri_downsampled',GS_fmri') ;
    %Fisher Z transformaiton
    GSCORR_fmri_normal(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr));
end

    parfor i = 1:nsims
        
        % GSCORR in firing rate
        GS_x_normal = mean(x,1) ;
        res_corr = corr(x',GS_x_normal') ;
        GSCORR_x_normal(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); %Fisher Z transformaiton
    
    
        % GSCORR in simulated bold
        sim_fmri_downsampled = downsample(sim_fmri',200)' ;
       % sim_fmri_hrf_downsampled = downsample(sim_fmri_hrf',200)' ;

        GS_fmri_normal = mean(sim_fmri_downsampled,1) ;
        res_corr = corr(sim_fmri_downsampled',GS_fmri_normal') ;
        GSCORR_fmri_normal(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); %Fisher Z transformaiton

       % GS_fmri_hrf_normal = mean(sim_fmri_hrf_downsampled,1) ;
       % res_corr = corr(sim_fmri_hrf_downsampled',GS_fmri_hrf_normal') ;
       % GSCORR_fmri_hrf_normal(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); %Fisher Z transformaiton
    
        % ACW
        for t=1:360 % Regions
       % [ACW0, ~, ~, ~] = ACW_kaan(sim_fmri_downsampled(t,:),1/dt,10,50);
       % [ACW0, ~, ~, ~] = ACW_kaan(sim_fmri_downsampled(t,:),1/2,10,50);

       % ACW with firing rate
       [ACW0, ~, ~, ~] = acw(x(t,:),100,false);
       ACW_x_normal(t,i) = ACW0 ;

       % Bw simulation
      % [ACW0, ~, ~, ~] = acw(sim_fmri_downsampled(t,:),1/2,false);
      % ACW_bw_normal(t,i) = ACW0 ;
       
       % HRF simulation
       % [ACW0, ~, ~, ~] = acw(sim_fmri_hrf_downsampled(t,:),1/2,false);
       % ACW_hrf_normal(t,i) = ACW0 ;

        end
        toc
        fprintf("\ni = %d\n", i)
    end

%% Figure for model
[x,time] = wilsoncowan_RK2(tau, b, W, k, s, C, tspan) ;
sim_fmri = bw(x,200,1/100) ;
sim_fmri_downsampled = downsample(sim_fmri',200)' ; % Downsampled at 0.5 Hz

f = figure ;
    set(gcf,'color','w');
    set(gcf,'units','normalized','position',[0 0 0.3 0.3])

subplot(1,3,1)
imdata = imread('C:/Users/kaan/Downloads/Wilcox.jpeg') ;    
imshow(imdata)
axis square
title('Wilson – Cowan Network','FontSize',20)

subplot(1,3,2)
plot(time,x')
axis square  
title('Simulated time-series')
set(gca, 'FontSize', 15);  % Adjust the font size as needed
ylabel('Firing rate', 'FontSize', 20);
xlim([100, 300])
xlabel('time(secs)')


subplot(1,3,3)
plot(sim_fmri_downsampled')
axis square  
set(gca, 'FontSize', 15);  % Adjust the font size as needed
title({'Balloon-Windkessel','Hemodynamic Modelling from firing rates'})
subtitle('Downsampled 0.5 Hz')
ylabel('BOLD simulation', 'FontSize', 20);
xlim([0, 91])

print('E:/EIB/FIGURES/FIG_4_A','-dpng','-r300');
%% Swarmchart for ACW and GSCORR
W = load('G:Drive''ım\Research\Myelin_Self\averageConnectivity_Fpt.mat'); %DTI common matrix
W = 10.^(W.Fpt); % Matrix Log10 olceginde. Kullanmadan once 10.^ yapmak gerekiyor
n = length(W);
W(logical(eye(n))) = 1; % Makes diagonal elemenst equal to 1. 
W(~logical(eye(n))) = (2*360)*(W(~logical(eye(n))) ./ sum(W(~logical(eye(n))), 'all'));
W_1 = log10(W) ; % Just to show a better plot

subplot(2,4,6)

imagesc(W_1)
axis square
colormap('jet');  % Set the colormap (you can choose a different one)
xlabel('ROI');
ylabel('');
yticks([]);
title("DTI Matrix (W)", 'FontSize',14)

% Get the diagonal elements
% diagonal_elements = W_1(logical(eye(n)));

axScale = axes('Position', [0.33, 0.1, 0.02, 0.35]); % Adjust the position as needed
imagesc(W_1(logical(eye(n))), 'Parent', axScale);
colormap('jet');
title('W_{ii} = 1','Interpreter', 'tex','FontSize',10)
xlabel('Recurrence','FontSize', 10)
ylabel('sum = 360','FontSize',10)
xticks([]);
yticks([]);

clim(axScale, [min(W_1(:)); max(W_1(:))] );

subplot(2,4,7)
% ACW firing rate
%{
% This is not true, compares NOISE because compares 30 simulation. Should
% be comparing the regions.
group1 = nanmean(ACW_x_normal(logical(DtiMyelin(:,1)-1),:),1);
group2 = nanmean(ACW_x_normal(~logical(DtiMyelin(:,1)-1),:),1);

labels = {'Self', 'Non-Self'};
data_labels = [repmat(labels(1), 30, 1); repmat(labels(2), 30, 1)];

plot_dat = table([group1,group2]', data_labels, 'VariableNames', {'values','class'}) ;

x = categorical(plot_dat.class, labels, 'Ordinal',true) ;
y = plot_dat.values ;
c = [ repmat(hex2rgb("E64B35"),30,1); repmat(hex2rgb("4DBBD5"),30,1) ] ;

swarmchart(x,y, 100,c, 'filled', 'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);

% Paired t-test
[h, ~, ~, stats] = ttest(group1, group2);

if h
    text(1.5, max([mean(group1), mean(group2)]) + 0.1, '*','FontSize', 20, 'HorizontalAlignment', 'center');
else
     text(1.5, max([mean(group1), mean(group2)]) + 0.1, 'ns','FontSize', 20, 'HorizontalAlignment', 'center');
end

set(gca, 'FontSize', 20);  % Adjust the font size as needed
ylabel('ACW', 'FontSize', 20);
title('Firing Rate','FontSize',25)
%}

group1 = nanmean(ACW_x_normal(logical(DtiMyelin(:,1)-1),:),2);
group2 = nanmean(ACW_x_normal(~logical(DtiMyelin(:,1)-1),:),2);

labels = {'Self', 'Non-Self'};
data_labels = [repmat(labels(1), length(group1), 1); repmat(labels(2), length(group2), 1)];

plot_dat = table([group1;group2], data_labels, 'VariableNames', {'values','class'}) ;

x = categorical(plot_dat.class,labels, 'Ordinal',true) ;
y = plot_dat.values ;
c = [ repmat(hex2rgb("E64B35"), length(group1),1); repmat(hex2rgb("4DBBD5"),length(group2),1) ] ;

swarmchart(x,log(y), 100,c, 'filled', 'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
hold on 
boxchart(x,log(plot_dat.values) ) ;

% t test
[h, p, ~, stats] = ttest2(group1, group2) ;

if h
    text(1.5, max(log([mean(group1), mean(group2)])) + 0.2, '*','FontSize', 30, 'HorizontalAlignment', 'center');
else
     text(1.5, max(log([mean(group1), mean(group2)])) + 0.2, 'ns','FontSize', 20, 'HorizontalAlignment', 'center');
end

% Effect-size 
mean1 = mean(group1);
mean2 = mean(group2);

std1 = std(group1);
std2 = std(group2);

n1 = length(group1);
n2 = length(group2);

pooledStd = sqrt(((n1 - 1) * std1^2 + (n2 - 1) * std2^2) / (n1 + n2 - 2));
effectSize = (mean1 - mean2) / pooledStd;

disp(['Cohen''s d: ', num2str(effectSize)]);

set(gca, 'FontSize', 15);  % Adjust the font size as needed
ylabel('log(ACW)', 'FontSize', 15);
title('Firing Rate','FontSize',20)
ylim([-0.55, 0.1])
subtitle(['Cohen''s d: ', num2str(effectSize)],"FontSize",15)


% GSCORR with Balloon
subplot(2,4,8)

group1 = nanmean( GSCORR_fmri_normal(logical(DtiMyelin(:,1)-1),:),2);
group2 = nanmean( GSCORR_fmri_normal(~logical(DtiMyelin(:,1)-1),:),2);

group1 = nanmean( GSCORR_fmri_normal(logical(DtiMyelin(:,1)-1),:),1)';
group2 = nanmean( GSCORR_fmri_normal(~logical(DtiMyelin(:,1)-1),:),1)';

labels = {'Self', 'Non-Self'};
data_labels = [repmat(labels(1), length(group1), 1); repmat(labels(2), length(group2), 1)];

plot_dat = table([group1;group2], data_labels, 'VariableNames', {'values','class'}) ;

x = categorical(plot_dat.class,labels, 'Ordinal',true) ;
y = plot_dat.values ;
c = [ repmat(hex2rgb("E64B35"), length(group1),1); repmat(hex2rgb("4DBBD5"),length(group2),1) ] ;

swarmchart(x,log(y), 100,c, 'filled', 'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
hold on 
boxchart(x,log(plot_dat.values) ) ;

% t test
[h, p, ~, stats] = ttest2(group1, group2) ;

if h
    text(1.5, max(log([mean(group1), mean(group2)])) + 0.6, '*','FontSize', 30, 'HorizontalAlignment', 'center');
else
     text(1.5, max(log([mean(group1), mean(group2)])) + 0.6, 'ns','FontSize', 20, 'HorizontalAlignment', 'center');
end

% Effect-size 
mean1 = mean(group1);
mean2 = mean(group2);

std1 = std(group1);
std2 = std(group2);

n1 = length(group1);
n2 = length(group2);

pooledStd = sqrt(((n1 - 1) * std1^2 + (n2 - 1) * std2^2) / (n1 + n2 - 2));
effectSize = (mean1 - mean2) / pooledStd;

disp(['Cohen''s d: ', num2str(effectSize)]);

set(gca, 'FontSize', 15);  % Adjust the font size as needed
ylabel('log(GSCORR)', 'FontSize', 15);
title(' Balloon-Windkessel model','FontSize',15)
ylim([-4.5, -1.5])
subtitle(['Cohen''s d: ', num2str(effectSize)],"FontSize",15)

%% Simulation with different diagonal 
% Assuming your first array is arr1 with length 360
arr1 = cell2mat(myDataParcels(:,3)) ;
% Calculate the inversely related values within the range of 0.3 to 2.5
arr2_inverse = 2 ./ arr1;

% Scale the values to be within the specified range with a minimum of 0.5
arr2_scaled = 0.3 + (arr2_inverse - min(arr2_inverse)) / (max(arr2_inverse) - min(arr2_inverse)) * (2.5 - 0.3);
arr1_scaled = 0.3 + (arr1 - min(arr1)) / (max(arr1) - min(arr1)) * (2.5 - 0.3);

% Display or use the new array arr2_scaled
disp(arr2_scaled); 

scaled_array_final = arr2_scaled * 360/sum(arr2_scaled) ; 
%%
%W(logical(eye(n))) = resultMatrix ;
%W(logical(eye(n))) = arr2_inverse ;
%W(logical(eye(n))) = arr2_scaled ;
W(logical(eye(n))) = scaled_array_final ;

W(~logical(eye(n))) = (2*360)*(W(~logical(eye(n))) ./ sum(W(~logical(eye(n))), 'all'));

%% Plot for myelin and recurrent connection 
subplot(2,4,5)
one_line = ones( 1, 360 ) ;
scatter(cell2mat(myDataParcels(cell2mat(myDataParcels(:,2))==1,3)),...  % myelin
     one_line(cell2mat(myDataParcels(:,2))==1),...
     30, ...
     hex2rgb("4DBBD5"),...
     "filled" )
axis square
hold on
scatter(cell2mat(myDataParcels(cell2mat(myDataParcels(:,2))==2,3)),...  % myelin
     one_line(cell2mat(myDataParcels(:,2))==2),...
     30, ...
     hex2rgb("E64B35"),...
     "filled" )
axis square

xlabel('Intracranial Myelin Content','FontSize',14);
ylabel('Recurrent Connections','FontSize',14);
title('No Scaling', 'Interpreter','tex','FontSize',14)

subplot(2,4,1)
scatter(cell2mat(myDataParcels(cell2mat(myDataParcels(:,2))==1,3)),...  % myelin
     scaled_array_final(cell2mat(myDataParcels(:,2))==1),...
     30, ...
     hex2rgb("4DBBD5"),...
     "filled" )
axis square
hold on
scatter(cell2mat(myDataParcels(cell2mat(myDataParcels(:,2))==2,3)),...  % myelin
     scaled_array_final(cell2mat(myDataParcels(:,2))==2),...
     30, ...
     hex2rgb("E64B35"),...
     "filled" )
axis square

xlabel('Intracranial Myelin Content','FontSize',14);
ylabel('Recurrent Connections','FontSize',14);
title('Inverse Scaling', 'Interpreter','tex','FontSize',14)

hold off


%% 
nsims = 100 ;
GSCORR_x_InvScaled= zeros(360,nsims);
GSCORR_fmri_InvScaled = zeros(360,nsims);
GSCORR_fmri_hrf_InvScaled = zeros(360,nsims);

ACW_x_InvScaled = zeros(360,nsims);
ACW_bw_InvScaled = zeros(360,nsims);
ACW_hrf_InvScaled = zeros(360,nsims);
    parfor i = 1:nsims
        tic
        %% Set up W
        %% run simulations
        [x,time] = wilsoncowan_RK2(tau, b, W, k, s, C, tspan);
        sim_fmri = bw(x,200,1/100) ;
        %sim_fmri_hrf = hrfconv(x, time, 10/1000)

        % GSCORR in firing rate
        GS_x = mean(x,1) ;
        res_corr = corr(x',GS_x') ;
        GSCORR_x_InvScaled(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); %Fisher Z transformaiton
    
    
        % GSCORR in simulated bold
        sim_fmri_downsampled = downsample(sim_fmri',200)' ;
        %sim_fmri_hrf_downsampled = downsample(sim_fmri_hrf',200)' ;

        GS_fmri = mean(sim_fmri_downsampled,1) ;
        res_corr = corr(sim_fmri_downsampled',GS_fmri') ;
        GSCORR_fmri_InvScaled(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); %Fisher Z transformaiton

        %GS_fmri_hrf = mean(sim_fmri_hrf_downsampled,1) ;
        %res_corr = corr(sim_fmri_hrf_downsampled',GS_fmri_hrf') ;
        %GSCORR_fmri_hrf_InvScaled(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); %Fisher Z transformaiton
    
        % ACW
        for t=1:360 % Regions
       % [ACW0, ~, ~, ~] = ACW_kaan(sim_fmri_downsampled(t,:),1/dt,10,50);
       % [ACW0, ~, ~, ~] = ACW_kaan(sim_fmri_downsampled(t,:),1/2,10,50);

       % ACW with firing rate
       [ACW0, ~, ~, ~] = acw(x(t,:),100,false);
       ACW_x_InvScaled(t,i) = ACW0 ;

       % Bw simulation
       % [ACW0, ~, ~, ~] = acw(sim_fmri_downsampled(t,:),1/2,false);
       % ACW_bw_InvScaled(t,i) = ACW0 ;
       
       % HRF simulation
       % [ACW0, ~, ~, ~] = acw(sim_fmri_hrf_downsampled(t,:),1/2,false);
       % ACW_hrf_InvScaled(t,i) = ACW0 ;

        end
        toc
        fprintf("\ni = %d\n", i)
    end
    %% Swarmchart for ACW and GSCORR
% Myelin is inversely scaled according to myelin part2, run after above
    subplot(2,4,2)
W(logical(eye(n))) = scaled_array_final ;
W(~logical(eye(n))) = (2*360)*(W(~logical(eye(n))) ./ sum(W(~logical(eye(n))), 'all'));
W_2 = log10(W) ; 

imagesc(W_2)
axis square
colormap('jet');  % Set the colormap (you can choose a different one)
xlabel('ROI');
ylabel('');
yticks([]);
title('DTI Matrix (W)', 'Interpreter','tex','FontSize',14)

axScale = axes('Position', [0.33, 0.58, 0.02, 0.35]); % Adjust the position as needed
imagesc(W_2(logical(eye(n))), 'Parent', axScale);
title('W_{ii} is inversely \newline scaled to myelin','FontSize',10, 'Interpreter','tex')
xlabel('Recurrence','FontSize', 10)
ylabel('sum = 360','FontSize',10)
xticks([]);
yticks([]);

clim(axScale, [min(W_2(:)); max(W_2(:))] );

    subplot(2,4,3)
% ACW firing rate
% Across all regions of simulations
% RUN FOR WITHOUT MEANING REGIONS IN SELF AND NONSELF
%SELF = ACW_x_InvScaled(logical(DtiMyelin(:,1)-1),:) ;
%NONSELF = ACW_x_InvScaled(~logical(DtiMyelin(:,1)-1),:) ;
%data_labels = [repmat(labels(1), 30, 1); repmat(labels(2), 30, 1)];

% % Number of colors
% numColors = 30; % Different colors for different simulations

% Generate random RGB values
% randomColors = rand(numColors, 3);

% Repeat each row 33 times and stack vertically
% repeatedColors_Self = kron(randomColors, ones(33, 1));
% repeatedColors_NonSelf = kron(randomColors, ones(327, 1));

% labels = {'Self', 'Non-Self'};
% data_labels = [repmat(labels(1), length(SELF(:)), 1); repmat(labels(2), length(NONSELF(:)), 1)];

% plot_dat = table([SELF(:);NONSELF(:)], data_labels, 'VariableNames', {'values','class'}) ;

% x = categorical(plot_dat.class,labels, 'Ordinal',true) ;
% y = plot_dat.values ;
% c = [ repmat(hex2rgb("E64B35"), length(SELF(:)),1); repmat(hex2rgb("4DBBD5"),length(NONSELF(:)),1) ] ;
% c = [ repeatedColors_Self ; repeatedColors_NonSelf ] ;

% swarmchart(x,y, 100,c, 'filled', 'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
%{
mean1 = mean(SELF(:));
mean2 = mean(NONSELF(:));

std1 = std(SELF(:));
std2 = std(NONSELF(:));

% Calculate Cohen's d
n1 = length(SELF(:));
n2 = length(NONSELF(:));

pooledStd = sqrt(((n1 - 1) * std1^2 + (n2 - 1) * std2^2) / (n1 + n2 - 2));
effectSize = (mean1 - mean2) / pooledStd;

disp(['Cohen''s d: ', num2str(effectSize)]);

%}
%  ACW firing rate
% Self and Nonself regions are averaged than compared. Each dot is a region
group1 = nanmean(ACW_x_InvScaled(logical(DtiMyelin(:,1)-1),:),2);
group2 = nanmean(ACW_x_InvScaled(~logical(DtiMyelin(:,1)-1),:),2);

labels = {'Self', 'Non-Self'};
data_labels = [repmat(labels(1), length(group1), 1); repmat(labels(2), length(group2), 1)];

plot_dat = table([group1;group2], data_labels, 'VariableNames', {'values','class'}) ;

x = categorical(plot_dat.class,labels, 'Ordinal',true) ;
y = plot_dat.values ;
c = [ repmat(hex2rgb("E64B35"), length(group1),1); repmat(hex2rgb("4DBBD5"),length(group2),1) ] ;

swarmchart(x,log(y), 100,c, 'filled', 'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
hold on 
boxchart(x,log(plot_dat.values) ) ;
% t test
[h, p, ~, stats] = ttest2(group1, group2) ;

if h
    text(1.5, max(log([mean(group1), mean(group2)])) + 0.2, '*','FontSize', 30, 'HorizontalAlignment', 'center');
else
     text(1.5, max(log([mean(group1), mean(group2)])) + 0.2, 'ns','FontSize', 20, 'HorizontalAlignment', 'center');
end

% Effect-size 
mean1 = mean(group1);
mean2 = mean(group2);

std1 = std(group1);
std2 = std(group2);

n1 = length(group1);
n2 = length(group2);

pooledStd = sqrt(((n1 - 1) * std1^2 + (n2 - 1) * std2^2) / (n1 + n2 - 2));
effectSize = (mean1 - mean2) / pooledStd;

disp(['Cohen''s d: ', num2str(effectSize)]);

set(gca, 'FontSize', 15);  % Adjust the font size as needed
ylabel('log(ACW)', 'FontSize', 15);
title('Firing Rate','FontSize',20)
ylim([-0.55, 0.55])
subtitle(['Cohen''s d: ', num2str(effectSize)],"FontSize",15)



subplot(2,4,4)

% GSCORR with Balloon

group1 = nanmean(GSCORR_fmri_InvScaled(logical(DtiMyelin(:,1)-1),:),2);
group2 = nanmean(GSCORR_fmri_InvScaled(~logical(DtiMyelin(:,1)-1),:),2);

labels = {'Self', 'Non-Self'};
data_labels = [repmat(labels(1), length(group1), 1); repmat(labels(2), length(group2), 1)];

plot_dat = table([group1;group2], data_labels, 'VariableNames', {'values','class'}) ;

x = categorical(plot_dat.class,labels, 'Ordinal',true) ;
y = plot_dat.values ;
c = [ repmat(hex2rgb("E64B35"), length(group1),1); repmat(hex2rgb("4DBBD5"),length(group2),1) ] ;

swarmchart(x, log(y), 100,c, 'filled', 'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
hold on 
boxchart(x,log(plot_dat.values) ) ;
% t test
[h, p, ~, stats] = ttest2(group1, group2) ;

if h
    text(1.5, max(log([mean(group1), mean(group2)]) ) + 0.2, '*','FontSize', 30, 'HorizontalAlignment', 'center');
else
     text(1.5, max(log([mean(group1), mean(group2)]) ) + 0.2, 'ns','FontSize', 20, 'HorizontalAlignment', 'center');
end

% Effect-size 
mean1 = mean(group1);
mean2 = mean(group2);

std1 = std(group1);
std2 = std(group2);

n1 = length(group1);
n2 = length(group2);

pooledStd = sqrt(((n1 - 1) * std1^2 + (n2 - 1) * std2^2) / (n1 + n2 - 2));
effectSize = (mean1 - mean2) / pooledStd;

disp(['Cohen''s d: ', num2str(effectSize)]);

set(gca, 'FontSize', 15);  % Adjust the font size as needed
ylabel('log(GSCORR)', 'FontSize', 15);
title(' Balloon-Windkessel model','FontSize',15)
ylim([-4.5,-1.5])
subtitle(['Cohen''s d: ', num2str(effectSize)],"FontSize",15)


print('E:/EIB/FIGURES/FIG_4','-dpng','-r300');
    
    %%

    GSCORR_x_InvScaled;
    GSCORR_fmri_InvScaled;
    GSCORR_fmri_hrf_InvScaled;

    ACW_x_InvScaled ;
    ACW_bw_InvScaled ;
    ACW_hrf_InvScaled ;

[H,P] = ttest2(GSCORR_x_InvScaled(~logical(DtiMyelin(:,1)-1),:), ...
    GSCORR_x_InvScaled(logical(DtiMyelin(:,1)-1),:)); % For each subject making a comparison. No results

[H,P] = ttest2(ACW_bw_InvScaled(~logical(DtiMyelin(:,1)-1),:), ...
    ACW_bw_InvScaled(logical(DtiMyelin(:,1)-1),:)); % For each subject making a comparison. No results

%% Null scenario
W(logical(eye(n))) = zeros(1,360) ;
W(~logical(eye(n))) = (2*360)*(W(~logical(eye(n))) ./ sum(W(~logical(eye(n))), 'all'));

nsims = 100 ;
GSCORR_x_ZeroW= zeros(360,nsims);
GSCORR_fmri_ZeroW = zeros(360,nsims);
GSCORR_fmri_hrf_ZeroW = zeros(360,nsims);

ACW_x_ZeroW = zeros(360,nsims);
ACW_bw_ZeroW = zeros(360,nsims);
ACW_hrf_ZeroW = zeros(360,nsims);
    parfor i = 1:nsims
        tic
        %% Set up W
        %% run simulations
        [x,time] = wilsoncowan_RK2(tau, b, W, k, s, C, tspan);
        sim_fmri = bw(x,200,1/100) ;
        %sim_fmri_hrf = hrfconv(x, time, 10/1000)

        % GSCORR in firing rate
        GS_x = mean(x,1) ;
        res_corr = corr(x',GS_x') ;
        GSCORR_x_ZeroW(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); %Fisher Z transformaiton
    
    
        % GSCORR in simulated bold
        sim_fmri_downsampled = downsample(sim_fmri',200)' ;
        %sim_fmri_hrf_downsampled = downsample(sim_fmri_hrf',200)' ;

        GS_fmri = mean(sim_fmri_downsampled,1) ;
        res_corr = corr(sim_fmri_downsampled',GS_fmri') ;
        GSCORR_fmri_ZeroW(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); %Fisher Z transformaiton

        %GS_fmri_hrf = mean(sim_fmri_hrf_downsampled,1) ;
        %res_corr = corr(sim_fmri_hrf_downsampled',GS_fmri_hrf') ;
        %GSCORR_fmri_hrf_InvScaled(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); %Fisher Z transformaiton
    
        % ACW
        for t=1:360 % Regions
       % [ACW0, ~, ~, ~] = ACW_kaan(sim_fmri_downsampled(t,:),1/dt,10,50);
       % [ACW0, ~, ~, ~] = ACW_kaan(sim_fmri_downsampled(t,:),1/2,10,50);

       % ACW with firing rate
       [ACW0, ~, ~, ~] = acw(x(t,:),100,false);
       ACW_x_ZeroW(t,i) = ACW0 ;

       % Bw simulation
       % [ACW0, ~, ~, ~] = acw(sim_fmri_downsampled(t,:),1/2,false);
       % ACW_bw_InvScaled(t,i) = ACW0 ;
       
       % HRF simulation
       % [ACW0, ~, ~, ~] = acw(sim_fmri_hrf_downsampled(t,:),1/2,false);
       % ACW_hrf_InvScaled(t,i) = ACW0 ;

        end
        toc
        fprintf("\ni = %d\n", i)
    end

    % Self and Nonself regions are averaged than compared. Each dot is a region
group1 = nanmean(ACW_x_ZeroW(logical(DtiMyelin(:,1)-1),:),2);
group2 = nanmean(ACW_x_ZeroW(~logical(DtiMyelin(:,1)-1),:),2);

labels = {'Self', 'Non-Self'};
data_labels = [repmat(labels(1), length(group1), 1); repmat(labels(2), length(group2), 1)];

plot_dat = table([group1;group2], data_labels, 'VariableNames', {'values','class'}) ;

x = categorical(plot_dat.class,labels, 'Ordinal',true) ;
y = plot_dat.values ;
c = [ repmat(hex2rgb("E64B35"), length(group1),1); repmat(hex2rgb("4DBBD5"),length(group2),1) ] ;

swarmchart(x,log(y), 100,c, 'filled', 'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
hold on 
boxchart(x,log(plot_dat.values) ) ;
% t test
[h, p, ~, stats] = ttest2(group1, group2) ;

if h
    text(1.5, max(log([mean(group1), mean(group2)])) + 0.2, '*','FontSize', 30, 'HorizontalAlignment', 'center');
else
     text(1.5, max(log([mean(group1), mean(group2)])) + 0.2, 'ns','FontSize', 20, 'HorizontalAlignment', 'center');
end

% Effect-size 
mean1 = mean(group1);
mean2 = mean(group2);

std1 = std(group1);
std2 = std(group2);

n1 = length(group1);
n2 = length(group2);

pooledStd = sqrt(((n1 - 1) * std1^2 + (n2 - 1) * std2^2) / (n1 + n2 - 2));
effectSize = (mean1 - mean2) / pooledStd;

disp(['Cohen''s d: ', num2str(effectSize)]);

set(gca, 'FontSize', 15);  % Adjust the font size as needed
ylabel('log(ACW)', 'FontSize', 15);
title('Firing Rate','FontSize',20)
ylim([-0.55, 0.55])
subtitle(['Cohen''s d: ', num2str(effectSize)],"FontSize",15)



subplot(2,4,4)

% GSCORR with Balloon

group1 = nanmean(GSCORR_fmri_ZeroW(logical(DtiMyelin(:,1)-1),:),2);
group2 = nanmean(GSCORR_fmri_ZeroW(~logical(DtiMyelin(:,1)-1),:),2);

labels = {'Self', 'Non-Self'};
data_labels = [repmat(labels(1), length(group1), 1); repmat(labels(2), length(group2), 1)];

plot_dat = table([group1;group2], data_labels, 'VariableNames', {'values','class'}) ;

x = categorical(plot_dat.class,labels, 'Ordinal',true) ;
y = plot_dat.values ;
c = [ repmat(hex2rgb("E64B35"), length(group1),1); repmat(hex2rgb("4DBBD5"),length(group2),1) ] ;

swarmchart(x, log(y), 100,c, 'filled', 'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
hold on 
boxchart(x,log(plot_dat.values) ) ;
% t test
[h, p, ~, stats] = ttest2(group1, group2) ;

if h
    text(1.5, max(log([mean(group1), mean(group2)]) ) + 0.2, '*','FontSize', 30, 'HorizontalAlignment', 'center');
else
     text(1.5, max(log([mean(group1), mean(group2)]) ) + 0.2, 'ns','FontSize', 20, 'HorizontalAlignment', 'center');
end

% Effect-size 
mean1 = mean(group1);
mean2 = mean(group2);

std1 = std(group1);
std2 = std(group2);

n1 = length(group1);
n2 = length(group2);

pooledStd = sqrt(((n1 - 1) * std1^2 + (n2 - 1) * std2^2) / (n1 + n2 - 2));
effectSize = (mean1 - mean2) / pooledStd;

disp(['Cohen''s d: ', num2str(effectSize)]);

set(gca, 'FontSize', 15);  % Adjust the font size as needed
ylabel('log(GSCORR)', 'FontSize', 15);
title(' Balloon-Windkessel model','FontSize',15)
ylim([-4.5,-1.5])
subtitle(['Cohen''s d: ', num2str(effectSize)],"FontSize",15)

