%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%

function merged_vp = merge_vp(vp1,vp2,args,im_args)

%initial new vp, just average
merged_vp.vp = [];
%initial new lsegs group
merged_vp.lsegs = [vp1.lsegs;vp2.lsegs];

merged_vp.var = 1e6;
merged_vp.tr_pt = 1e6;
merged_vp.pan_var = 1e6;
merged_vp.tilt_var = 1e6;
merged_vp.consistency_measure = [];
merged_vp.gauss_error = [];
merged_vp.mean_gauss_error = -1;
merged_vp.pan_tilt = [];
merged_vp.gauss_point = [];



%filter out outlier lsegs
% p = p_of_lseg_given_v_probability(merged_vp.lsegs,merged_vp,args.stdev);
% [max_p,index]=max(p,[],2);
% chosen_lines_index = find(max_p > args.min_consistency / 2);
% merged_vp.lsegs = merged_vp.lsegs(chosen_lines_index,:); 

%% refine new vp
merged_vp = optimize_vp(merged_vp, args, [0,0], im_args.focal);


