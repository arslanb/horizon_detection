%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%

function lsegs = pre_process_lsegs(lsegs)

[n_lsegs, temp] = size(lsegs);
if temp < 4
    lsegs = [];
    return;
end

x1 = [lsegs(:, [1 3]) ones(size(lsegs, 1), 1)];
x2 = [lsegs(:, [2 4]) ones(size(lsegs, 1), 1)];

% angel
lsegs(:,5)= atan2( lsegs(:,3) - lsegs(:,4), lsegs(:,1)-lsegs(:,2) );

lsegs(:,7) = sqrt(sum((x1-x2).^2,2));

% line equations
l = cross(x1, x2);
l = l ./ repmat(sqrt(sum(l(:,1:2).^2,2)), 1, 3);
lsegs(:,8:10) = l; % nx3
pp = [0,0,1]';
r = l*pp; % nx1
lsegs(:,6) = r; % r
lsegs(:,11)=[1:n_lsegs]'; % ID