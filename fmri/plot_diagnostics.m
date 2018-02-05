
% Show the figure?
showfig = 'on';
save = 1;

% Read general diagnostics files:
censor  = csvread('censor_combined.1D');
outlier = csvread('outcount.1D')*100;
motion  = csvread('motion_enorm.1D');

% Want to plot only censored data, and at y=0.
censor(censor == 1) = NaN;

% X-axis of the plots
TRs = 1:length(censor);

% Regions of interest
rois = {'R_IFG','R_preSMA','R_ACC','R_GPe','R_GPi','R_STN',...
'L_ACC','R_Thal','R_NAcc','R_Caudate'};

% Average timeseries and fit plots for each ROI.
for roi = rois
    roi = roi{1};

    % Some files contain no data because 0 voxels survive
    % the F-statistic mask. These return read errors.
    try
        ave = csvread(['',roi,'_ave.1D']);
        fit = csvread(['',roi,'_fitts.1D']);
    catch error
        disp(['No data for ROI: ',roi])
        continue
    end

    % Don't pop-up on cluster.
    figure('visible', showfig)

    % Average and fit plot.
    hold on
    plot(TRs, [ave,fit])
    plot(TRs, censor, 'og')
    hold off
    ylabel('z-score')
    xlabel('TR [2.2s]')
    title([roi, ' Average BOLD Signal'], 'Interpreter', 'None')
    legend({'Observed','Fit','Censor'})
    info = ['R^2 = ', num2str(corr(ave,fit))];
    text(50,-1,info);
    grid('on')
    ylim([-2,2])

    % Make it bigger and save.
    set(gcf, 'Position', [100, 208, 1352, 740])
    if save; saveas(gcf,['plots/fits_', roi, '.png']); end
end


% Motion, censoring, and outliers plot.
figure('visible', showfig)

hold on
plot(TRs, censor, 'og')
h = plotyy(TRs, outlier, TRs, motion);
ylabel(h(1), 'Outliers [% of Voxels]')
ylabel(h(2), 'Motion Enorm [mm]')
xlabel('TR [2.2s]')
title('Quality')
grid('on')
hold off
legend({'Censor','Outlier','Motion'})
set(h(1), 'YLim', [0,6])
set(h(1), 'YTick', [0,1,2,3,4,5,6])
set(h(2), 'YLim', [0,3])
set(h(2), 'YTick', [0,0.5,1,1.5,2,2.5,3])

% Make it bigger and save.
% set(gcf, 'Position', [100, 208, 1352, 740])
if save; saveas(gcf,['plots/motion_', roi, '.png']); end

% Regressors
nregressors = 14;
regressors = NaN(TRs(end), nregressors);

% I would love to read X.nocensor.mat.1D but that's a huge pain.
regressors(:,1) = csvread('ideals/ideal_go_success.1D'    );
regressors(:,2) = csvread('ideals/ideal_stop_success.1D'  );
regressors(:,3) = csvread('ideals/ideal_stop_failure.1D'  );
regressors(:,4) = csvread('ideals/ideal_go_failure.1D'    );
regressors(:,5) = csvread('ideals/ideal_go_too_late.1D'   );
regressors(:,6) = csvread('ideals/ideal_go_wrong_key.1D'  );
regressors(:,7) = csvread('ideals/ideal_stop_too_early.1D');

% Add motion breakdown
motion_params = dlmread('motion_demean.1D', ' ');
motion_params = motion_params(:,2:end);
regressors(:,8:13) = motion_params;

% Add total motion
regressors(:,14) = motion;

% Name them
regnames = {'Go Success', 'Stop Success', 'Stop Failure', 'Go Failure', 'Go Too Late', ...
    'Go Wrong Key', 'Stop Too Early', 'Roll', 'Pitch', 'Yaw', 'dS', 'dL', 'dP', 'Motion'};

%
figure('visible', showfig)

for regnum = 1:7
    subplot(7,2,regnum*2 -1)
    plot(regressors(:,regnum))
    title(regnames{regnum})
end

for index = 1:7
    regnum = 7+index;
    subplot(7,2, index*2)
    plot(regressors(:,regnum))
    title(regnames{regnum})
end

% Make it bigger and save.
set(gcf, 'Position', [95, 1, 1396, 960])
if save; saveas(gcf,['plots/regressors.png']); end

