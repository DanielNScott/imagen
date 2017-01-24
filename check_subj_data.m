% Check SST Data

ss_presented(ss_presented == -999) = NaN;
ssd(ssd == -999) = NaN;
rt(rt == -999) = NaN;
inhibited(inhibited == -999) = NaN;

max_len = 0;
for i = 1:subj_idx(end)
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

for i = 100:100 %subj_idx(end)
   subj_msk = subj_idx == i;
   
   if sum(subj_msk) == 0;
      continue
   end
   
   subj_ss  = ss_presented(subj_msk);
   subj_ssd = ssd(subj_msk);
   subj_rt  = rt(subj_msk);
   subj_inh = inhibited(subj_msk);

   figure();
   subplot(3,1,1)
   plot(subj_ss, '.')
   title(['subj ' num2str(i) ' ss'])
   
   subplot(3,1,2)
   plot(subj_ssd, '.')
   title(['subj ' num2str(i) ' ssd'])
   
   subplot(3,1,3)
   plotyy(1:length(subj_rt), subj_rt, 1:length(subj_rt), subj_inh)
   title(['subj ' num2str(i) ' rt and inhib'])
   
   %subplot(2,2,4)
   %plot(subj_inh, '.')
   %title(['subj ' num2str(i) ' inh'])
end