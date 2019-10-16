%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%

function direction = check_direction(lseg, vp)

e1 = lseg([1,3]);
e2 = lseg([2,4]);

direction = norm(e1-vp(1:2)) > norm(e2-vp(1:2));