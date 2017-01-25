% This file checks the SST data found in a .csv file with the
% format specified in the BEESTS program.
% 
% To use it, load such a csv file into the appropriate 5 vectors:
%  - inhibited
%  - rt
%  - ss_presented
%  - ssd
%  - subj_idx

plot_flag = 1;
suspect_ids = [175, 183];

% Convert missing data to NaNs:
rt(rt == -999) = NaN;
ssd(ssd == -999) = NaN;
inhibited(inhibited == -999) = NaN;
ss_presented(ss_presented == -999) = NaN;

max_len = 0;
for i = suspect_ids(1):1:suspect_ids(2)
   subj_msk = subj_idx == i;
   subj_rt = rt(subj_msk);
   subj_ss = ss_presented(subj_msk);
   
   max_len = max(length(subj_rt),max_len);
   
   if length(subj_rt) > 600;
      disp(i)
   end
   
   
   n_non_ss_rts = length(subj_ss(subj_ss == 0));
   
   if n_non_ss_rts > 0 && n_non_ss_rts <= 2;
      disp(['Too few non SS RTs for subject id: ' num2str(i)])
   end
end

for i = suspect_ids(1):suspect_ids(2)
   subj_msk = subj_idx == i;
   
   if sum(subj_msk) == 0;
      continue
   end
   
   subj_ss  = ss_presented(subj_msk);
   subj_ssd = ssd(subj_msk);
   subj_rt  = rt(subj_msk);
   subj_inh = inhibited(subj_msk);

   tf = sum(~subj_ss & isnan(subj_rt));
   disp(['Potential Observed TFs for Subject ', num2str(i), ': ', num2str(tf)])
   disp(['Number of Observed RTs for Subject ', num2str(i), ': ', num2str(length(~isnan(subj_rt)))])
   
   if plot_flag
      figure();
      subplot(3,2,1)
      plot(subj_ss, '.')
      title(['subj ' num2str(i) ' ss'])

      subplot(3,2,3)
      plot(subj_ssd, '.')
      title(['subj ' num2str(i) ' ssd'])

      subplot(3,2,5)
      plotyy(1:length(subj_rt), subj_rt, 1:length(subj_rt), subj_inh)
      title(['subj ' num2str(i) ' rt and inhib'])

      edges = 0:50:1000;
      
      subplot(3,2,2)
      N = histc(subj_rt(subj_inh == 0), edges);
      bar(edges,N,'histc')
      title(['subj ' num2str(i) ' successful inhib histogram'])
      
      subplot(3,2,4)
      N = histc(subj_ssd, edges);
      bar(edges,N,'histc')
      title(['subj ' num2str(i) ' ssd histogram'])

      subplot(3,2,6)
      N = histc(subj_rt, edges);
      bar(edges,N,'histc')
      title(['subj ' num2str(i) ' rt histogram'])
   end
end