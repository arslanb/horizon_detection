%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%

function vp = min_error_vp(vp,args,pp,focal)
% Given vp cluster, principle point and focal length, compute the least
% error vp;
% Input: vp - vanishing point and clustering strcture
%        pp - priciple point
%        focal - focal length    

stdev=args.stdev;
%% if vp.vp is empty, then use least square to estimate it.
if isempty(vp.vp)
    vp.vp = least_square_vp(vp.lsegs);
end

%% use max-log to refine
vp = max_log_vp(vp,args,pp,focal);
consistency_measure = p_of_lseg_given_v_probability(vp.lsegs,vp.vp,args.stdev);
[p,sorted_index]=sort(consistency_measure,'descend');
if consistency_measure(args.min_num_lsegs) < max(consistency_measure)*0.5
    back_bone_index = sorted_index(1:args.min_num_lsegs);
else
    back_bone_index = find(consistency_measure > max(consistency_measure)*0.5);% & consistency_measure > args.min_consistency);
end
temp_lsegs=vp.lsegs;
vp.lsegs=vp.lsegs(back_bone_index,:);

%% construct trace vector
n = size(vp.lsegs,1); % number of lsegs;
w=zeros(n*(n-1)/2,1);  % weight vector;
tr=zeros(n*(n-1)/2,1);  % trace;
pan_var = zeros(n*(n-1)/2,1);
tilt_var = zeros(n*(n-1)/2,1);
cross_vps = zeros(n*(n-1)/2,3);


i=0;
for j = 1:n-1
    for k = j+1:n
        cross_vp = cross(vp.lsegs(j,8:10), vp.lsegs(k,8:10));
%         if cross_vp(3)==0
%             i = i+1;
%             cross_vp(3) = 1e-6;
%             cross_vp = cross_vp / cross_vp(3);
%             cross_vps(i,:) = cross_vp;
%             tr(i) = 1e6;
%            pan_var(i) = 0.5e6;
%            tilt_var(i) = 0.5e6;
%         end
        cross_vp = cross_vp / cross_vp(3);
        if check_direction(vp.lsegs(j,1:4),cross_vp) == check_direction(vp.lsegs(j,1:4),vp.vp)  && ... 
           check_direction(vp.lsegs(k,1:4),cross_vp) == check_direction(vp.lsegs(k,1:4),vp.vp) 
           
           J_img2pt = jacobian_image_to_pt(cross_vp, focal);
           J_line2img = jacobian_line_pair_to_image( vp.lsegs(j,1:4), vp.lsegs(k,1:4) );
           J = J_img2pt*J_line2img;
           cov_pt = J*J';
           t = trace(cov_pt); % can ignore the component caused by endpoint noise, assuming they are all the same.
           if isnan(t)
               continue;
           end
           i = i+1;
           cross_vps(i,:) = cross_vp;
           tr(i) = t;
           pan_var(i) = cov_pt(1,1);
           tilt_var(i) = cov_pt(2,2);
        end
    end
end

if i<length(w)
    w(i+1:end)=[];
    tr(i+1:end)=[];
    cross_vps(i+1:end,:) = [];
    pan_var(i+1:end) = [];
    tilt_var(i+1:end) = [];
end


%% compute optimal weights
recip_tr = ones(size(tr))./tr;
denominator = sum(recip_tr);
w = recip_tr ./ denominator;

%% compute least-error vp
ww = w.^2;
vp.vp = w' * cross_vps;
vp.vp = vp.vp / vp.vp(3);
vp.tr_pt = ww'*tr;
vp.pan_var = ww'*pan_var;
vp.tilt_var = ww'*tilt_var;
vp.pan_tilt = img2pt(vp.vp,pp,focal);
vp.gauss_point = img2gaussPoint(vp.vp,pp,focal);
%vp.lsegs=temp_lsegs;



