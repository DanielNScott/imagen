% check_param_consistency

path = '/home/dan/projects/imagen/data/model fits/';
set1 = readtable([path, fname1], 'Delimiter', ';');
set2 = readtable([path, fname2], 'Delimiter', ';');

params = {'mu_go_subj_', 'mu_stop_subj_', 'tau_go_subj_', 'tau_stop_subj_', ...
     'sigma_go_subj_', 'sigma_stop_subj_'}; 

n_subj = 198;
n_params = length(params);
n_diagnost = 8;
diagnostics = NaN(n_subj, n_params*n_diagnost);
for subj = 1:198
    for param_num = 1:n_params
        name = [params{param_num}, num2str(subj)];
        
        try
            avg1 = mean(set1.(name));
            avg2 = mean(set2.(name));
            
            std1 = std(set1.(name));
            std2 = std(set2.(name));

            diff_avg = abs(avg1 - avg2) / ((abs(avg1) + abs(avg2))/2);
            diff_std = abs(std1 - std2) / ((abs(std1) + abs(std2))/2);

            % The computation of these gets completed after all subjs are added. 
            diff_avg_rel = abs(avg1 - avg2);
            diff_std_rel = abs(std1 - std2);

            tmp = [avg1, avg2, std1, std2, diff_avg, diff_std, diff_avg_rel, diff_std_rel];
            
            slice_inds = ((param_num-1)*n_diagnost + 1):(param_num*n_diagnost);
            diagnostics(subj, slice_inds) = tmp;
            
        catch message
            continue
            
        end
    end
end

% 
avg_diff_inds = 7:n_diagnost:n_diagnost*n_params;
avg1_inds = 1:n_diagnost:n_diagnost*n_params;
avg2_inds = 2:n_diagnost:n_diagnost*n_params;

mad_avg1 = mad(diagnostics(:, avg1_inds), 1);
mad_avg2 = mad(diagnostics(:, avg2_inds), 1);

diff_avg_rel = 100* diagnostics(:, avg_diff_inds) ./ ((mad_avg1 + mad_avg2)/2);
diagnostics(:, avg_diff_inds) = diff_avg_rel;

%
std_diff_inds = 8:n_diagnost:n_diagnost*n_params;
std1_inds = 3:n_diagnost:n_diagnost*n_params;
std2_inds = 4:n_diagnost:n_diagnost*n_params;

mad_std1 = mad(diagnostics(:, std1_inds), 1);
mad_std2 = mad(diagnostics(:, std2_inds), 1);

diff_std_rel = 100* diagnostics(:, std_diff_inds) ./ ((mad_std1 + mad_std2)/2);
diagnostics(:, std_diff_inds) = diff_std_rel;

% Normalize average differences relative to inter-individual variation

figure()
set(gcf, 'Position', [100, 435, 1309, 529])

subplot(2,2,1)
inds = 5:n_diagnost:n_diagnost*n_params;
semilogy(1:198, diagnostics(:, inds), 'o');
legend(params, 'Location', 'SouthEast', 'Interpreter', 'None')
ylabel('difference [%]')
xlabel('subject number')
title('Chain Mean Differences')
set(gca,'YGrid', 'on')

subplot(2,2,2)
inds = 6:n_diagnost:n_diagnost*n_params;
semilogy(1:198, diagnostics(:, inds), 'o');
legend(params, 'Location', 'SouthEast', 'Interpreter', 'None')
ylabel('difference [%]')
xlabel('subject number')
title('Chain Std. Dev. Differences')
set(gca,'YGrid', 'on')

subplot(2,2,3)
inds = 7:n_diagnost:n_diagnost*n_params;
semilogy(1:198, diagnostics(:, inds), 'o');
legend(params, 'Location', 'SouthEast', 'Interpreter', 'None')
ylabel('difference as % MAD [%]')
xlabel('subject number')
title({'Chain Mean Differences', 'Relative to Between Subj. Variation'})
set(gca,'YGrid', 'on')

subplot(2,2,4)
inds = 8:n_diagnost:n_diagnost*n_params;
semilogy(1:198, diagnostics(:, inds), 'o');
legend(params, 'Location', 'SouthEast', 'Interpreter', 'None')
ylabel('difference as % MAD [%]')
xlabel('subject number')
title({'Chain Std. Dev. Differences ', 'Relative to Between Subj. Variation'})
set(gca,'YGrid', 'on')

