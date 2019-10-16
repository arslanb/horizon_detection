%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%

function lsegs = sort_lsegs(lsegs,vp)

v_to_1 = lsegs(:,[1,3])-repmat(vp(1:2),size(lsegs,1),1);
v_to_2 = lsegs(:,[2,4])-repmat(vp(1:2),size(lsegs,1),1);

need_reverse_index = find(sum(v_to_1.*v_to_1,2) < sum(v_to_2.*v_to_2,2));
temp = lsegs(need_reverse_index,[1,3]);
lsegs(need_reverse_index,[1,3])= lsegs(need_reverse_index,[2,4]);
lsegs(need_reverse_index,[2,4])= temp;
% one_to_two = lsegs(:,[1,3])-lsegs(:,[2,4]);
% cos_angle_1 = sum(v_to_1 .* one_to_two);
% cos_angle_1 = sum(v_to_2 .* one_to_two);
       