%plot_betas

betas = roibetas14(:,3:end);
names = betas.Properties.VariableNames;

for i = 1:length(names)
   
   subplot(5,10,i)
   histogram( table2array(betas(:,i)) )
   title(names(i))
    
end