%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%
% return cost of a vp with respect to a set of line segments sorted_lsegs
% assume sorted_lsegs(i,:) = [x1,x2,y1,y2] and [x1,y1] is furhter from vp
% than [x2,y2] and all line segments are completely on one side of the vp,
% i.e., 
function cost = vpCostFun(vp,sorted_lsegs)

n_lines = size(sorted_lsegs,1);

x1=sorted_lsegs(:,1); x2=sorted_lsegs(:,2); y1=sorted_lsegs(:,3); y2=sorted_lsegs(:,4);
x=repmat(vp(1),n_lines,1); y=repmat(vp(2),n_lines,1);
L=sorted_lsegs(:,7);

temp = (y1-y2).*x + (x2-x1).*y + x1.*y2 - x2.*y1;

u=sqrt( (x-x1).^2+(y-y1).^2-temp.^2 ) ./ L;
v=temp./L;
   
cost = sum(L.^2.*v.^2 ./ (u.^2 + (u-L).^2));




