%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%
% Detect horizon of a given image
% Input:    im_fn: input image name (with full directory)
%           options:  algorithm parameters and user options. Leaving it
%           empty then default values will be used.

% Output:   h: [h1,h2,h3], a 3x1 array, the horizon line represented in the homogeneous
%           coordinate system. x*h1+y*h2+h3 = 0.  (0,0) is the center of the image.
%           It will return h as an empty array when the algorithm does not find a horizon.
% 
% Authr: Yiliang Xu, Kitware

function h = horizon_detection(im_fn, options)

close all;

%% parameters
if ~(exist('args','var') && ~isempty(options) )
    options.debug = 0;
    options.opt_option = 'min-error';  % choose algorithm variations for vanishing point detection: 'min-error' or 'max-log' 
    options.color = [0,0,0;
        0,0,1;
        0,1,0;
        1,0,0;
        1,1,0;
        1,0,1;
        0,1,1;
        0.5, 0.5, 0;
        0.5, 0, 0.5;
        0, 0.5, 0.5;
        0, 0, 0.5;
        1, 0, 0.5;
        ];
    options.sigma_gate = 1; % how many std away from mean counts as outlier, the smaller the stricter.
    options.min_num_lsegs = 6;  % minimum size of line segment group
    options.min_intersect = 1; % min intersection between two lsegs groups that is considered to be merged.
    options.stdev = 1.5; % line segment endpoint noise standard deviation.
    options.min_consistency = 1/(sqrt(2*pi)*options.stdev) * exp(-options.sigma_gate^2/2);
    options.lsegs_generator = 'lsd'; % method for line segment detection.
    options.lsd_param = 0.3; % parameter for LSD
    options.min_lseg_length = 0.03; % minimum accepted length of line segment w.r.t. to the image diagnal length
    options.min_agular_difference = 10; % minimum difference in orientation to distinguish two different line segment groups 
    %options.min_vp_prob_diff = 0.2; 
    options.nIteration = 10; % number of EM iterations
end


%% read image
im = imread(im_fn);
[H,W,c] = size(im);

%% image dependent parameters
if ~(exist('im_args','var') && ~isempty(im_args) )
    im_args.focal_range = [0.15, 5] * max(H,W);
    im_args.focal = max(H,W);
end

% minimum length of line segments
min_len = norm([W,H])*options.min_lseg_length;

%% detect line segments
if strcmp(char(options.lsegs_generator), 'lsd') % extract long strgight line segments using lsd c++ exe
    lsegs = extract_lsd(im, min_len, options.lsd_param);
else
    lsegs = extract_lsd(im, min_len, options.lsd_param);
end

num_lsegs = size(lsegs,1);
if size(lsegs,1) < options.min_num_lsegs
    h=[];
    disp(['Find only ',num2str(num_lsegs), ' line segments, less than required ', num2str(options.min_num_lsegs), '. Return empty horizon line']);
    return;
end

% visualize all line segments
if options.debug
    figure;
    imshow(im);
    hold on;
    for j = 1:size(lsegs,1)
        plot(lsegs(j,1:2)+W/2, lsegs(j,3:4)+H/2,'-r');
        hold on;
    end
    title('all detected line segments');
    hold off;
end

% change line segments coordinate and create other fields
lsegs = pre_process_lsegs(lsegs);


%% detect vanishing points
[vps, prob] = vp_probability_EM(lsegs, im, options, im_args);
if isempty(prob)
    h=[];
    disp('J-linkage did not find any hypothesized vanishing point, return empty horizon line');
    return;
end
fh_grouping=show_vps_and_lsegs(vps,im,[]);


%% generate horizon
[h,zenith,v_vps,h_vps] = horizon(vps,[H,W],[0,0],im_args.focal_range,options.opt_option);


if ~isempty(h)
    % optional: to visualize horizon in image
    [fh_horizon] = draw_lines_in_img(im,h);
else
    disp('Did not find any horizon');
end
    

disp('Done horizon detection');









