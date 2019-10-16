%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%

function vp = max_log_vp(vp,args,pp,focal)

stdev=args.stdev;
%% if vp.vp is empty, then use least square to estimate it.
if isempty(vp.vp)
    vp.vp = least_square_vp(vp.lsegs);
end

sorted_lsegs = sort_lsegs(vp.lsegs, vp.vp(1:2));

[x,fval] = fminsearch(@(x) vpCostFun(x,sorted_lsegs), vp.vp(1:2));
v = [x,1];
vp.vp=v;

vp.consistency_measure = p_of_lseg_given_v_probability(vp.lsegs,vp.vp,stdev); % not most efficient
% remove = find(vp.consistency_measure < args.min_consistency);
% vp.lsegs(remove,:)=[];
% vp.consistency_measure(remove)=[];
vp.gauss_error = -log(vp.consistency_measure * sqrt(2*pi)*stdev) / (2*stdev^2); % not most efficient
vp.mean_gauss_error = mean(vp.gauss_error);
vp.pan_tilt = img2pt(vp.vp,pp,focal);
vp.gauss_point = img2gaussPoint(vp.vp,pp,focal);
