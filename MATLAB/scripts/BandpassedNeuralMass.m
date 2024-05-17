source_path = '/home/kaansocat/MYELIN' ;
addpath(strcat(source_path,'/MATLAB/func/') );
W = load(strcat(source_path,'/DATA/averageConnectivity_Fpt.mat') ); %DTI common matrix
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
load(strcat(source_path,'/DATA/myDataParcels.mat'),'myDataParcels') % For later

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

[bp, ap] = cheby1(N, Rp, Wp, 'bandpass');

nsims = 1000;

GSCORR_fmri_normal = zeros(360,nsims);
GSCORR_fmri_normal_bp = zeros(360,nsims);
GSCORR_fmri_normal_dt = zeros(360,nsims);
GSCORR_fmri_normal_bp_dt = zeros(360,nsims);

ACW_x_normal = zeros(360,nsims);
%ACW_x_normal_bp = zeros(360,nsims);
%ACW_x_normal_bp_dt = zeros(360,nsims);

GSCORR_fmri_InvScaled = zeros(360,nsims);
GSCORR_fmri_InvScaled_bp = zeros(360,nsims);
GSCORR_fmri_InvScaled_dt = zeros(360,nsims);
GSCORR_fmri_InvScaled_bp_dt = zeros(360,nsims);

ACW_x_InvScaled = zeros(360,nsims);
%ACW_x_InvScaled_bp = zeros(360,nsims);
%ACW_x_InvScaled_bp_dt = zeros(360,nsims);
%% Simulation with Wii = 1 
W = load(strcat(source_path,'/DATA/averageConnectivity_Fpt.mat') ); %DTI common matrix
W = 10.^(W.Fpt); % Matrix Log10 olceginde. Kullanmadan once 10.^ yapmak gerekiyor
n = length(W);
W(logical(eye(n))) = 1; % Makes diagonal elemenst equal to 1. 
W(~logical(eye(n))) = (2*360)*(W(~logical(eye(n))) ./ sum(W(~logical(eye(n))), 'all'));

%% run simulations
% Those simulations saves the each wilson-cowan simulations.
% However, size is too big. Thus measures are calculated wihtin each loop
% in below section


%sim_fmri_subj = zeros(360, 18002, nsims);
%sim_rate_subj = zeros(360, 20002, nsims);

%{
parfor i = 1:nsims
        tic
        [x,time] = wilsoncowan_RK2(tau, b, W, k, s, C, tspan);
        sim_fmri = bw(x,200,1/100) ;
       % sim_fmri_hrf = hrfconv(x, time, 10/1000)
        sim_rate_subj(:,:,i) = x ;
        sim_fmri_subj(:,:,i) = sim_fmri ;
        toc
        fprintf("\ni = %d\n", i)
end

% ACW0
parfor i=1:nsims
    for t=1:360
       % ACW with firing rate
       [ACW0, ~, ~, ~] = acw(sim_rate_subj(t,:,i),100,false);
       ACW_x_InvScaled(t,i) = ACW0 ;
    end
end

% Only Bandpassed 
parfor i=1:nsims
    % Downsampling fmri simulation to match with HCP data. TR = 2 sampling
    % is sampling the 0.5 Hz. Thus our TR = 0.72, thus we will sampling at
    % Fs (1/Tr). dt = 100 in simulation, so .72 x 100 = 72 should be
    % downsampling rate.
    sim_fmri_downsampled = downsample(sim_fmri_subj(:,:,i)', 72)' ;
    % Bandpassing the data
    filteredData = zeros(size(sim_fmri_downsampled));
    % Apply the filter to each ROI (each voxel/ROI)
    for t = 1:size(sim_fmri_downsampled, 1)
        filteredData(t, :) = filtfilt(bp, ap, sim_fmri_downsampled(t, :));
    end
    
    % Calculating GS from filtered Data
    GS_fmri = mean(sim_fmri_downsampled,1) ;
    GS_fmri_filtered = filtfilt(bp, ap, GS_fmri);
    res_corr = corr(filteredData',GS_fmri_filtered') ;
    %Fisher Z transformaiton
    GSCORR_fmri_InvScaled(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); 
end

% Bandpassed and detrended
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
        filteredData(t, :) = filtfilt(bp, ap, sim_fmri_downsampled_dt(t, :));
    end
    
    % Calculating GS from filtered Data
    GS_fmri = mean(sim_fmri_downsampled_dt,1) ;
    GS_fmri_filtered = filtfilt(bp, ap, GS_fmri);
    res_corr = corr(filteredData',GS_fmri_filtered') ;
    %Fisher Z transformaiton
    GSCORR_fmri_InvScaled(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); 
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
    GSCORR_fmri_InvScaled(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr));
end

%}
        %% run simulations
parfor i = 1:nsims
    tic
    [x,time] = wilsoncowan_RK2(tau, b, W, k, s, C, tspan);
    sim_fmri = bw(x,200,1/100) ;

    % Downsampling fmri simulation to match with HCP data. TR = 2 sampling
    % is sampling the 0.5 Hz. Thus our TR = 0.72, thus we will sampling at
    % Fs (1/Tr). dt = 100 in simulation, so .72 x 100 = 72 should be
    % downsampling rate.

    %sim_fmri_downsampled = downsample(sim_fmri', 72)' ;
    sim_fmri_downsampled = downsample(sim_fmri', 200)' ;

    %=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    % Directly saving, without any other processing 
    %=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    % Calculating GS 
    GS_fmri_normal = mean(sim_fmri_downsampled,1) ;
    res_corr = corr(sim_fmri_downsampled',GS_fmri_normal') ;
    %Fisher Z transformaiton
    GSCORR_fmri_normal(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); 

    %=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    %       O N L Y   B A N D P A S S I N G
    %=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    filteredData = zeros(size(sim_fmri_downsampled));
    sim_fmri_downsampled_trans = sim_fmri_downsampled' ;
    % Apply the filter to each ROI (each voxel/ROI)
    for t = 1:size(sim_fmri_downsampled, 1)
        filteredData(t, :) = filtfilt(bp, ap, sim_fmri_downsampled_trans(:, t));
    end
    
    % Calculating GS from filtered Data
    GS_fmri = mean(sim_fmri_downsampled,1) ;
    GS_fmri_filtered = filtfilt(bp, ap, GS_fmri) ;
    res_corr = corr(filteredData',GS_fmri_filtered') ;
    % Fisher Z transformaiton
    GSCORR_fmri_normal_bp(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); 

    %=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    %       O N L Y   D E T R E N D I N G
    %=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    sim_fmri_downsampled_detrended  = detrend(sim_fmri_downsampled')' ;
    % Calculating GS 
    GS_fmri_normal = mean(sim_fmri_downsampled_detrended,1) ;
    res_corr = corr(sim_fmri_downsampled_detrended',GS_fmri_normal') ;
    %Fisher Z transformaiton
    GSCORR_fmri_normal_dt(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); 
    
    %=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    % D E T R E N D I N G  +  B A N D P A S S I N G
    %=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    sim_fmri_downsampled_detrended_trans = sim_fmri_downsampled_detrended' ;
    filteredData = zeros(size(sim_fmri_downsampled));
   
    for t = 1:size(sim_fmri_downsampled, 1)
        filteredData(t, :) = filtfilt(bp, ap, sim_fmri_downsampled_detrended_trans(:, t)) ;
    end
    
    % Calculating GS from filtered Data
    GS_fmri = mean(sim_fmri_downsampled_detrended,1) ;
    GS_fmri_filtered = filtfilt(bp, ap, GS_fmri);
    res_corr = corr(filteredData',GS_fmri_filtered') ;
    % Fisher Z transformaiton
    GSCORR_fmri_normal_bp_dt(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); 


    % ACW
    for t=1:360 % Regions
   % ACW with firing rate
   [ACW0, ~, ~, ~] = acw(x(t,:),100,false);
   ACW_x_normal(t,i) = ACW0 ;
    end
    toc
    fprintf("\ni = %d\n", i)
end

%{
figure;
hold on;
histogram(GSCORR_fmri_normal(:,i), 'Normalization', 'pdf', 'FaceColor', [0.255, 0.412, 0.882], 'FaceAlpha', 0.5);
histogram(GSCORR_fmri_normal_dt(:,i), 'Normalization', 'pdf', 'FaceColor', [0.863, 0.078, 0.235], 'FaceAlpha', 0.5);
histogram(GSCORR_fmri_normal_bp(:,i), 'Normalization', 'pdf', 'FaceColor', [0.133, 0.545, 0.133], 'FaceAlpha', 0.5);
histogram(GSCORR_fmri_normal_bp_dt(:,i), 'Normalization', 'pdf', 'FaceColor',[0.855, 0.647, 0.125], 'FaceAlpha', 0.5);
hold off;
%}
%% Simulation with different diagonal 
W = load(strcat(source_path,'/DATA/averageConnectivity_Fpt.mat') ); %DTI common matrix
W = 10.^(W.Fpt); % Matrix Log10 olceginde. Kullanmadan once 10.^ yapmak gerekiyor
n = length(W);
% Assuming your first array is arr1 with length 360
arr1 = cell2mat(myDataParcels(:,3)) ;
% Calculate the inversely related values within the range of 0.3 to 2.5
arr2_inverse = 2 ./ arr1;

% Scale the values to be within the specified range with a minimum of 0.5
arr2_scaled = 0.3 + (arr2_inverse - min(arr2_inverse)) / (max(arr2_inverse) - min(arr2_inverse)) * (2.5 - 0.3);

scaled_array_final = arr2_scaled * 360/sum(arr2_scaled) ; 

W(logical(eye(n))) = scaled_array_final ;

W(~logical(eye(n))) = (2*360)*(W(~logical(eye(n))) ./ sum(W(~logical(eye(n))), 'all'));

        %% run simulations
parfor i = 1:nsims
    tic
    [x,time] = wilsoncowan_RK2(tau, b, W, k, s, C, tspan);
    sim_fmri = bw(x,200,1/100) ;

    % Downsampling fmri simulation to match with HCP data. TR = 2 sampling
    % is sampling the 0.5 Hz. Thus our TR = 0.72, thus we will sampling at
    % Fs (1/Tr). dt = 100 in simulation, so .72 x 100 = 72 should be
    % downsampling rate.

    %sim_fmri_downsampled = downsample(sim_fmri', 72)' ;
    sim_fmri_downsampled = downsample(sim_fmri', 200)' ;

    %=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    % Directly saving, without any other processing 
    %=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    % Calculating GS 
    GS_fmri_normal = mean(sim_fmri_downsampled,1) ;
    res_corr = corr(sim_fmri_downsampled',GS_fmri_normal') ;
    %Fisher Z transformaiton
    GSCORR_fmri_InvScaled(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); 

    %=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    %       O N L Y   B A N D P A S S I N G
    %=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

    filteredData = zeros(size(sim_fmri_downsampled));
    sim_fmri_downsampled_trans = sim_fmri_downsampled' ;
    % Apply the filter to each ROI (each voxel/ROI)
    for t = 1:size(sim_fmri_downsampled, 1)
        filteredData(t, :) = filtfilt(bp, ap, sim_fmri_downsampled_trans(:, t));
    end
    
    % Calculating GS from filtered Data
    GS_fmri = mean(sim_fmri_downsampled,1) ;
    GS_fmri_filtered = filtfilt(bp, ap, GS_fmri) ;
    res_corr = corr(filteredData',GS_fmri_filtered') ;
    % Fisher Z transformaiton
    GSCORR_fmri_InvScaled_bp(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); 

    %=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    %       O N L Y   D E T R E N D I N G
    %=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    sim_fmri_downsampled_detrended  = detrend(sim_fmri_downsampled')' ;
    % Calculating GS 
    GS_fmri_normal = mean(sim_fmri_downsampled_detrended,1) ;
    res_corr = corr(sim_fmri_downsampled_detrended',GS_fmri_normal') ;
    %Fisher Z transformaiton
    GSCORR_fmri_InvScaled_dt(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); 
    
    %=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    % D E T R E N D I N G  +  B A N D P A S S I N G
    %=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
    sim_fmri_downsampled_detrended_trans = sim_fmri_downsampled_detrended' ;
    filteredData = zeros(size(sim_fmri_downsampled));
   
    for t = 1:size(sim_fmri_downsampled, 1)
        filteredData(t, :) = filtfilt(bp, ap, sim_fmri_downsampled_detrended_trans(:, t)) ;
    end
    
    % Calculating GS from filtered Data
    GS_fmri = mean(sim_fmri_downsampled_detrended,1) ;
    GS_fmri_filtered = filtfilt(bp, ap, GS_fmri);
    res_corr = corr(filteredData',GS_fmri_filtered') ;
    % Fisher Z transformaiton
    GSCORR_fmri_InvScaled_bp_dt(:,i) =  0.5 * log((1 + res_corr) ./ (1 - res_corr)); 


    % ACW
    for t=1:360 % Regions
   % ACW with firing rate
   [ACW0, ~, ~, ~] = acw(x(t,:),100,false);
   ACW_kaan(x(t,:),100,100,false )
   ACW_x_InvScaled(t,i) = ACW0 ;
    end
    toc
    fprintf("\ni = %d\n", i)
end

script_name = mfilename('fullpath');
description = ['This data called as V2, because has a different downsampling' ...
    'to Tr = 2, not to TR = 0.71. This is to compare results with the latest' ...
    'version of the paper'];
disp('Finished');
save(strcat(source_path,'/DATA/1000SubjSim_v2.mat') )
%% Figure 
group1 = nanmean(GSCORR_fmri_normal_bp_dt(logical(DtiMyelin(:,1)-1),:),1)';
group2 = nanmean(GSCORR_fmri_normal_bp_dt(~logical(DtiMyelin(:,1)-1),:),1)';

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

% 

group1 = nanmean(ACW_x_InvScaled(logical(DtiMyelin(:,1)-1),:),1)';
group2 = nanmean(ACW_x_InvScaled(~logical(DtiMyelin(:,1)-1),:),1)';

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






