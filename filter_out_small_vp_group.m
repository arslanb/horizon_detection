%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%

function vps = filter_out_small_vp_group(vps,min_num_lsegs)

remove = [];
if min_num_lsegs < 2
   min_num_lsegs = 2; 
end


for i = 1:length(vps)
   if size(vps(i).lsegs, 1) < min_num_lsegs
       remove(end+1) = i;
       disp(['removing vp: ', num2str(i)]);
   end
end

vps(remove) = [];