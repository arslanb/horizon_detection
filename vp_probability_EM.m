%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%
% Compute vanishing points 
% Input:
% lsegs( [x1 x2 y1 y2 angle r length line_equation(3x1) ID] )
% im: image
% Output:
% vps(1xn_vp structs) - the found vanishing points; 
% p(nlsegs, n_vp) - the confidence that each line belongs to each vp
%
% Authr: Yiliang Xu, Kitware

function [vps, p] = vp_probability_EM(lsegs, im, options, im_args)


%% prepare
[H,W,c] = size(im);
imsize = [H,W];
nlsegs = size(lsegs, 1);



%% J-linkage clustering
vps = j_linkage_vp2(im,lsegs,options,[]);
fh=show_vps_and_lsegs(vps,im,options.color);
ngroups = length(vps);
if ngroups<1
    p=[];
    return;
end


%% After initialization, filter out small clusters
vps = filter_out_small_vp_group(vps,options.min_num_lsegs);
ngroups = length(vps);
if ngroups<1
    p=[];
    return;
end
if options.debug
    show_vps_and_lsegs(vps,im,options.color);
end


%% initialize EM
for i = 1:length(vps)
    lsegs_length = size(vps(i).lsegs,1)
    
    vps(i).vp = least_square_vp(vps(i).lsegs);
    vps(i).var = 1e6;
    vps(i).tr_pt = 1e6;
    vps(i).pan_var = 1e6;
    vps(i).tilt_var = 1e6;
    vps(i).consistency_measure = [];
    vps(i).gauss_error = [];
    vps(i).mean_gauss_error = -1;
    vps(i).pan_tilt = [];
    vps(i).gauss_point = [];
        
    % use cost function minimization to refine vp
    vps(i) = optimize_vp(vps(i),options, [0,0], im_args.focal);
    
end


%% do EM
p = p_of_lseg_given_v_probability(lsegs,vps,options.stdev);
for iter = 1:options.nIteration
    old_vps = vps;
    vps = []; 
    oldp = p;
    p = p_of_lseg_given_v_probability(lsegs,old_vps,options.stdev);    
    
    % for each vp
    new_vp_idx = 0;
    for i = 1:length(old_vps);
        % E step
        [max_p,index]=max(p,[],2);
        chosen_lines_index = find(index==i & max_p > options.min_consistency );
        if length(chosen_lines_index) < options.min_num_lsegs
            continue;
        end
        new_vp_idx = new_vp_idx + 1;
        vps(new_vp_idx,1).vp = old_vps(i).vp;
        vps(new_vp_idx,1).lsegs = lsegs(chosen_lines_index,:);
        vps(new_vp_idx,1).var = 1e6;
        vps(new_vp_idx,1).tr_pt = 1e6;
        vps(new_vp_idx,1).pan_var = 1e6;
        vps(new_vp_idx,1).tilt_var = 1e6;
        vps(new_vp_idx,1).consistency_measure = [];
        vps(new_vp_idx,1).gauss_error = [];
        vps(new_vp_idx,1).mean_gauss_error = -1;
        vps(new_vp_idx,1).pan_tilt = [];
        vps(new_vp_idx,1).gauss_point = [];
        
        % M-step
        % use cost function minimization to further refine vp
        vps(new_vp_idx) = optimize_vp(vps(new_vp_idx),options, [0,0],im_args.focal);
    end

    if options.debug
        show_vps_and_lsegs(vps,im,options.color);
    end
    
    
    %disp('start to split and merge');
    % split and merge vps
    p = p_of_lseg_given_v_probability(lsegs,vps,options.stdev);
    vps = merge_split_vps2(vps,p,lsegs, options, im_args);
    if options.debug
        show_vps_and_lsegs(vps,im,options.color);
    end
    
    % merge vps as a whole
    p = p_of_lseg_given_v_probability(lsegs,vps,options.stdev);
    [vps,p] = merge_vps_angle_based(vps,p,options,im_args);
    if options.debug
        show_vps_and_lsegs(vps,im,options.color);
    end
    % merge vps as a whole
    p = p_of_lseg_given_v_probability(lsegs,vps,options.stdev);
    [vps,p] = merge_vps_angle_based(vps,p,options,im_args);
    if options.debug
        show_vps_and_lsegs(vps,im,options.color);
    end
     % merge vps as a whole
    p = p_of_lseg_given_v_probability(lsegs,vps,options.stdev);
    [vps,p] = merge_vps_angle_based(vps,p,options,im_args);
    if options.debug
        show_vps_and_lsegs(vps,im,options.color);
    end
  
    
end
