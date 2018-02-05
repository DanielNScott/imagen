function [] = plot_beta_hists(statsTable, nrois)

% Convenience
betas = statsTable(:,3:end);

% Get the variable names
names = statsTable.Properties.VariableNames;
names = names(3:end);

% Need a big window
figure();
set(gcf, 'Position', [22, 85, 1807, 858]);

% Look at each variable
for i = 1:length(names)
    
   % Ignore columns with no values
   col = table2array(betas(:,i));
   if all(col == 'NA')
       continue
   end
   disp(' ')
   disp(['Output for ', names{i}])

   % Some ROIs have no significant voxels, even though
   % others for the same subject do. These cases produce
   % exact zeros in output.
   disp(['Num. zero data: ', num2str(sum(col == 0))])
   disp(['Num. NaN data: ', num2str(sum(isnan(col)))])
   col(col == 0) = NaN;
   
   % Median, Max Abs. Deviation, and outlier mask 
   med = nanmedian(col);
   md  = mad(col,1);
   msk = abs(col) > med + 2.5*md;
   
   % Robust distribution
   robust = col;
   robust(msk) = NaN;
   
   % Get and output statistics
   [h, p, ci] = ttest(robust);
   disp(['  H: ', num2str(h)])
   disp(['  p: ', num2str(p)])
   disp([' CI: ', num2str(ci(1)),' ', num2str(ci(2))])
   disp(['Avg: ', num2str(nanmean(robust)) ])
   disp(['Std: ', num2str(nanstd(robust )) ])
   
   % Plotting
   nrows = ceil(size(betas,2)/nrois);
   
   subplot(nrows, nrois, i)
   hold on
   histogram(robust)
   ylims = get(gca,'YLim');
   line([ci(1), ci(1)], [ylims(1), ylims(2)], 'Color', 'r', 'LineStyle', '--')
   line([ci(2), ci(2)], [ylims(1), ylims(2)], 'Color', 'r', 'LineStyle', '--')
   set(gca,'YLim',ylims)
   hold off
   
   % Title the top rows
   roi_name  = strsplit(names{i}, '_');
   condition = [roi_name{2:end}];
   roi_name  = roi_name{1};

   column_num  = mod(i-1, nrois) + 1;
   if i <= nrois
       title(roi_name, 'Interpreter', 'None')
   else
       % Make sure the column header is 'correct'
       header_name = strsplit(names{column_num}, '_');
       header_name = header_name{1};
       if header_name ~= roi_name
           error(['Column layout issue! ', i, roi_name, header_name])
       end
   end
   if column_num == 1
       ylabel(['\bf{',condition,'}'])
   end
   
end
end