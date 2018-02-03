
% Read general diagnostics files:
censor  = csvread('results/censor_s01_combined_2.1D');
outlier = csvread('results/outcount.r01.1D')*100;
motion  = csvread('results/motion_s01_enorm.1D');

% Want to plot only censored data, and at y=0.
censor(censor == 1) = NaN;

% X-axis of the plots
TRs = 1:length(censor);

% Regions of interest
rois = {'R_IFG','R_preSMA','R_AAL_ACC','R_GPe','R_GPi','R_STN',...
'L_AAL_ACC','R_Thalamus_AAL','R_NAcc','R_Caudate_AAL'};

% Average timeseries and fit plots for each ROI.
for roi = rois
    roi = roi{1};
    
    % Some files contain no data because 0 voxels survive 
    % the F-statistic mask. These return read errors.
    try
        ave = csvread(['results/',roi,'_ave.1D']);
        fit = csvread(['results/',roi,'_fitts.1D']);
    catch error
        disp(['No data for ROI: ',roi])
    end
    
    % Don't pop-up on cluster.
    figure('visible','off')
    
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
    saveas(gcf,['results/plots/fits_', roi, '.png'])
end


% Motion, censoring, and outliers plot.
figure('visible','off')

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
ylim([-2,2])

% Make it bigger and save.
set(gcf, 'Position', [100, 208, 1352, 740])
saveas(gcf,['results/plots/motion_', roi, '.png'])
