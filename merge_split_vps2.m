%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%

function vps = merge_split_vps2(vps,p,lsegs,args,im_args)

stdev = args.stdev;
min_consistency = args.min_consistency;
min_intersect = args.min_intersect;
min_num_lsegs = args.min_num_lsegs;

    
associate_matrix = p > min_consistency;

vps_assoc = cell(1,length(vps));
for i=1:length(vps)
    vps_assoc{i} = find(associate_matrix(:,i) == 1);
end


from=0; to=0; intersection = [];
max_intersect = min_intersect;
for i = 1:length(vps)-1
    for j = i+1:length(vps)
        inter = intersect(vps(i).lsegs(:,11), vps_assoc{j});
        overlap = length(inter); %/ length(vps(i).lsegs(:,11));
        if overlap > min_intersect
            if overlap > max_intersect
                from = i; to = j; max_intersect = overlap; intersection = inter;
            end
        end
        
    end
end


if from~=0 && to~=0 && max_intersect > min_intersect && ~isempty(intersection) 
    remaining_in_from = setdiff(vps(from).lsegs(:,11), intersection );
    remaining_in_to = setdiff(vps(to).lsegs(:,11), intersection );
    %length_remaining_in_from = length(remaining_in_from)
    %length_remaining_in_to = length(remaining_in_to)
    
    if length(remaining_in_to) < min_num_lsegs && length(remaining_in_from) < min_num_lsegs
        vps(to)=[];
        vps(from)=[];
    elseif length(remaining_in_to) < min_num_lsegs && length(remaining_in_from) >= min_num_lsegs
        vps(from).lsegs = lsegs(remaining_in_from,:);
        vps(from) = optimize_vp( vps(from), args, [0,0], im_args.focal );
        vps(to)=[];
    elseif length(remaining_in_to) >= min_num_lsegs && length(remaining_in_from) < min_num_lsegs
        vps(to).lsegs = lsegs(remaining_in_to,:);
        vps(to) = optimize_vp( vps(to), args, [0,0], im_args.focal );
        vps(from)=[];
    else
        vps(from).lsegs = lsegs(remaining_in_from,:);
        vps(from) = optimize_vp( vps(from), args, [0,0], im_args.focal );
        vps(to).lsegs = lsegs(remaining_in_to,:);
        vps(to) = optimize_vp( vps(to), args, [0,0], im_args.focal );
    end
    
end

%% if overlap part is big enough to be an independent group
if length(intersection) >= args.min_num_lsegs
    vps(end+1).lsegs = lsegs(intersection,:);
    vps(end).var = 1e6;
    vps(end).tr_pt = 1e6;
    vps(end).pan_var = 1e6;
    vps(end).tilt_var = 1e6;
    vps(end).consistency_measure = [];
    vps(end).gauss_error = [];
    vps(end).mean_gauss_error = -1;
    vps(end).pan_tilt = [];
    vps(end).gauss_point = [];
    vps(end) = optimize_vp( vps(end), args, [0,0], im_args.focal );
end

