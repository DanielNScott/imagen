
for i = 1:2;
   ss_presented = [ss_presented; ss_presented];
   inhibited    = [inhibited; inhibited];
   ssd          = [ssd; ssd];
   rt           = [rt; rt];
   
   subj_idx_tmp = subj_idx;
   
   for j = 1:subj_idx(end)
      subj_idx_tmp(subj_idx_tmp == j) = subj_idx(end) + j;
   end
      
   subj_idx = [subj_idx; subj_idx_tmp];
end

save_mat = [subj_idx, ss_presented, inhibited, ssd, rt];
csvwrite('fake_data_group.csv',save_mat)