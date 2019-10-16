%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%

function [vps,p] = merge_vps_angle_based(vps,p,args,im_args)

min_angular_difference = 10000; i_index = 0; j_index = 0;
for i = 1:length(vps)
    for j = i+1:length(vps)
        angular_difference = norm(vps(i).gauss_point-vps(j).gauss_point)/2;
        if angular_difference < min_angular_difference
            min_angular_difference = angular_difference; 
            i_index = i; 
            j_index = j;
        end
    end
end

if min_angular_difference < args.min_agular_difference / 180 * pi
    vps(i_index) = merge_vp(vps(i_index),vps(j_index), args, im_args);
    vps(j_index) = [];
    p(:, j_index) = [];
end