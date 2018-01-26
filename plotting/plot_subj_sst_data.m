function [] = plot_subj_sst_data(sst_data, subj_num, invisible)

subj = sst_data(sst_data(:,1) == subj_num,:);
subj(subj == -999) = NaN;

% Unpack data
ss_presented = logical(subj(:,2));
inhib = subj(:,3);
ssd   = subj(:,4);
rt    = subj(:,5);

% Number trials
ntrials = length(rt);
trial = 1:ntrials;


% Subplot controls
nrows      = 3;
ncols_bar  = 4;
ncols_time = 1;
nplts_time = 1;
plot_num   = 0;

% Plot all response times
figure()
if invisible
    set(gcf, 'Visible','off') % For use with looping & saving.
end
y_limits = [0,1000];

plot_num = plot_num + 1;
subplot(nrows, ncols_time, plot_num)
rt_plot = plot(trial, rt, 'o');

title('Responses')
xlabel('Trial Number')
ylabel('Reaction Time [ms]')
hold on

% Plot all errors of commission
eoc_inds = find(ss_presented & ~isnan(rt));
eoc_plot = plot(eoc_inds, rt(eoc_inds), 'or');

% Plot all omitted responses
no_response      = isnan(rt);
no_response_inds = find(no_response);

no_resp_line = plot(trial, ssd, 'om');
%Y = [100, 900];
%for i = 1:length(no_response_inds)
%    X = repmat(no_response_inds(i),2,1);
%    no_resp_line = line(X, Y);
%end

% Plot errors of omission
eoo_inds = find(no_response & ~ss_presented);

Y = [0, 1000];
eoo_line = line([NaN, NaN], [NaN, NaN], 'Color', [1,0,0,0.4]);
for i = 1:length(eoo_inds)
    X = repmat(eoo_inds(i),2,1);
    eoo_line = line(X, Y, 'Color', [1,0,0,0.4]);
end

npts = 7;
avg = movmean(rt, npts, 'omitnan');
avg_plt = plot(avg, 'Color', 'k');

ylim(y_limits)
xlim([trial(1), trial(end)])

set(gca, 'XGrid', 'on')
set(gca, 'YGrid', 'on')
hold off

legend([rt_plot, eoc_plot, no_resp_line, eoo_line, avg_plt] ...
    , {'Hit', 'Stop Fail', 'Stop Delay', 'Miss', 'Average'} ...
    , 'Location', 'NorthWest')

% Plot moving mean and moving SD
%plot_num = plot_num + 1;
%subplot(nrows, ncols_time, plot_num)
%npts = 7;
%avg = movmean(rt, npts, 'omitnan');
%std = movstd( rt, npts, 'omitnan');
%
%yu = (avg + 2*std)';
%yl = (avg - 2*std)';
%
%hold on
%avg_plt  = plot(avg, 'Color', 'm');
%std_fill = fill([trial, fliplr(trial)], [yu, fliplr(yl)], 'g', 'linestyle', 'none', 'FaceAlpha', 0.25);
%
%set(gca, 'XGrid', 'on')
%set(gca, 'YGrid', 'on')
%hold off
%
%title('Moving Averages')
%xlabel('Trial Number')
%ylabel('Reaction Time [ms]')
%
%legend([avg_plt, std_fill], {'7-pt Mean', '7-pt Std'}, 'Location', 'NorthWest')
%ylim(y_limits)
%xlim([trial(1), trial(end)])

% Histograms and bar plots:
plot_num = nplts_time*ncols_bar + 1;
subplot(nrows, ncols_bar, plot_num)
histogram(ssd(isnan(rt)))
title('Successful Stop SSDs')
xlabel('time [ms]')
ylabel('count')

plot_num = plot_num + 1;
subplot(nrows, ncols_bar, plot_num)
histogram(ssd(eoc_inds))
title('Error of Comission SSDs')
xlabel('time [ms]')
ylabel('count')


plot_num = plot_num + 1;
subplot(nrows, ncols_bar, plot_num)
histogram(rt(eoc_inds) - ssd(eoc_inds))
title('Error of Comission RT - SSD')
xlabel('time [ms]')
ylabel('count')

% Go trial reaction times
plot_num = plot_num + 1;
subplot(nrows, ncols_bar, plot_num)
histogram(rt(~ss_presented))
title('Go Trial RTs')
xlabel('time [ms]')
ylabel('count')

% Compute probability of failure to stop by SSD
bin_lower_lims = 0:50:500;

inhibs    = sum(sum((ssd(inhib == 1)) > bin_lower_lims, 2) == 1:11);
failures  = sum(sum((ssd(inhib == 0)) > bin_lower_lims, 2) == 1:11);
fail_prob = failures./ sum([inhibs; failures]);

plot_num = plot_num + 1;
subplot(nrows, ncols_bar, plot_num)

hold on
bar(bin_lower_lims, fail_prob, 'hist')
set(get(gca,'Children'),'FaceColor', [0, 0.4470, 0.7410])
set(get(gca,'Children'),'FaceAlpha', 0.6)

weights   = mnrfit(ssd(ss_presented), inhib(ss_presented)+1);
fit_pts   = 0:10:500;
fit_curve = glmval(weights, fit_pts, 'logit');
plot(fit_pts, fit_curve);
hold off

xlim([-50,550])
xticks([0:100:500])
ylim([0,1])
title('Empirical P(stop failure)')
xlabel('SSD [ms]')
ylabel('probability')

% Box-plots of stop-aligned RTs
maxlag = 7;
bin = zeros(ntrials,1);
for i = 2:ntrials
    if ~ss_presented(i)
        bin(i)  = bin(i-1) + 1;
    end
end

plot_num = plot_num + 1;
subplot(nrows, ncols_bar, plot_num)
binmatrix = repmat(rt,1,maxlag);
binmatrix(bin ~= 1:maxlag) = NaN;
boxplot(binmatrix, 'positions', 1:maxlag, 'labels', 1:maxlag)
set(gca, 'XGrid', 'on')
set(gca, 'YGrid', 'on')
xlim([0.5,maxlag+0.5])
ylim([100,900])
title('RT by time since SS')
xlabel('Consecutive go trials')
ylabel('RT [ms]')

hold on
glob_med = nanmedian(rt);
eoc_med  = nanmedian(rt(eoc_inds));
line([0, maxlag+1], [glob_med, glob_med], 'LineStyle', '--', 'Color', 'b')
line([0, maxlag+1], [eoc_med, eoc_med], 'LineStyle', '-.', 'Color', 'r')
hold off

% Increases and decreases
glob_std = nanstd(rt);
diff = nan( sum(bin == 1), maxlag-1);
for i = 2:maxlag
    msk = bin == i;
    npts = sum(msk);
    diff(1:npts, i-1)  = rt(msk) - rt( logical([msk(2:end);0]) );
end
diff = diff/glob_std;

plot_num = plot_num + 1;
subplot(nrows, ncols_bar, plot_num)
boxplot(diff, 'positions', 2:maxlag, 'labels', 2:maxlag)

set(gca, 'XGrid', 'on')
set(gca, 'YGrid', 'on')
title('\Delta RT relative to prev. trial')
xlabel('Consecutive go trials')
ylabel('coeff. of variation []')

% Calculate some SSRT measures, such as SSRT_med
% Ref:
% Band, G. P., Van Der Molen, M. W., & Logan, G. D. (2003).
% Horse-race model simulations of the stop-signal procedure.
% Acta Psychologica, 112(2), 105–142.
%
% Here I use logistic regression rather than linear regression to establish
% 50% response rate SSD.

plot_num = plot_num + 1;
subplot(nrows, ncols_bar, plot_num)
SSRT_med = nanmedian(rt(~ss_presented)) - ...
    mean(fit_pts((fit_curve - 0.5).^2 < 0.005));

bar(1, SSRT_med)
set(get(gca,'Children'),'FaceColor', [0, 0.4470, 0.7410])
set(get(gca,'Children'),'FaceAlpha', 0.6)
ylim([0,500])

% Widen figure
set(gcf,'Position', [1           1        1920         995])

end
