%plot_betas

betas = roibetas14q0(:,3:end);
names = betas.Properties.VariableNames;

for i = 1:length(names)
   
   if all(table2array(betas(1:5,i))== 'NA')
       continue
   end
   
   subplot(5,10,i)
   newbetas = table2array(betas(:,i));
   newmean = mean(newbetas);
   sd   = std(newbetas);
   
   %newbetas(abs(newbetas) > median(newbetas) + 2*mad(newbetas,1)) = NaN;
   
   histogram( newbetas )
   
   title(names(i), 'Interpreter', 'None')
   
  
   [h,p,ci] = ttest(table2array(betas(:,i)));
   disp(' ')
   disp(names{i})
   disp(['H=',num2str(h)])
   disp(['p=',num2str(p)])
   disp(['CIs', num2str(ci(1)),' ', num2str(ci(2))])
   disp(newmean)
   disp(sd)
end