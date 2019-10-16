%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%
% Extract line segments using LSD executable
% Input:    im: the input image.
%           min_len:  minimum length of a line segment.
%           lsd_param:  the -d parameter for LSD executable, which controls
%           the criteria for accepting a line segment. The higher, the more line segments will be returned.
% Output:   lsegs:  the Nx6 array representing all detected line segments.
%           Each row for one line segments: [x1 x2 y1 y2 orientation length]
%
% Aurthor: Yiliang Xu, Kitware

function lsegs = extract_lsd(im, min_len, lsd_param)


[H,W,c]=size(im);

imwrite(im,'input_pgm.pgm');
system_call = ['LSD_VC.exe -d ' num2str(lsd_param), ' input_pgm.pgm', ' ./lsd_output.txt'];
[status, result] = system(system_call)
lsegs = load('./lsd_output.txt');

%% LSD format
% % x1, y1, x2, y2, width, p, -log_nfa.
% For example, the line:
% 
%   159.232890 134.369601 160.325338 105.613616 2.735466 0.125000 17.212465
% 
% means that a line segment starting at point (159.232890,134.369601),
% ending at point (160.325338 105.613616) and of width 2.735466 was
% detected. An angle precision p of 0.125 was used, which means a
% gradient angle tolerance of p*180 = 0.125*180 = 22.5 degree. The
% opposite of the logarithm in base 10 of the NFA value of the detection
% was -log_10(NFA)=17.212465, so the NFA value was 10^(-17.2124656),
% roughly 6e-18. The length unit is the pixel and the origin of
% coordinates is the center of the top-left pixel (0,0).



%% change to appropriate format
% lines([x1 x2 y1 y2 angle r])

lsegs(:,6) = sqrt( (lsegs(:,1)-lsegs(:,3)).^2 + (lsegs(:,2)-lsegs(:,4)).^2 );
idx = find(lsegs(:,6) > min_len);
lsegs = lsegs(idx,:);

temp = lsegs(:,2);
lsegs(:,2) = lsegs(:,3);
lsegs(:,3) = temp;

N = size(lsegs,1);
lsegs(:,1:4) = lsegs(:,1:4) - repmat([W/2,W/2,H/2,H/2],N,1);
