%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%

function J = jacobian_line_pair_to_image( varargin )

if nargin == 2
   lseg1 = varargin{1};
   lseg2 = varargin{2};
   u_j1 = lseg1(1);
   v_j1 = lseg1(3);
   u_j2 = lseg1(2);
   v_j2 = lseg1(4);
   u_k1 = lseg2(1);
   v_k1 = lseg2(3);
   u_k2 = lseg2(2);
   v_k2 = lseg2(4);
elseif nargin == 4
   lseg1 = [varargin{1}(1),varargin{2}(1), varargin{1}(2),varargin{2}(2)];
   lseg2 = [varargin{3}(1),varargin{4}(1), varargin{3}(2),varargin{4}(2)];
   u_j1 = lseg1(1);
   v_j1 = lseg1(3);
   u_j2 = lseg1(2);
   v_j2 = lseg1(4);
   u_k1 = lseg2(1);
   v_k1 = lseg2(3);
   u_k2 = lseg2(2);
   v_k2 = lseg2(4); 
elseif nargin == 8
   u_j1 = varargin{1};
   v_j1 = varargin{3};
   u_j2 = varargin{2};
   v_j2 = varargin{4};
   u_k1 = varargin{1};
   v_k1 = varargin{3};
   u_k2 = varargin{2};
   v_k2 = varargin{4}; 
else
   error('Incorrect input for jacobian_line_pair_to_image'); 
end

J = [ ((u_k1 - u_k2)*(v_j1 - v_j2)*(u_j2*v_k1 - u_k1*v_j2 - u_j2*v_k2 + u_k2*v_j2 + u_k1*v_k2 - u_k2*v_k1))/(u_j1*v_k1 - u_k1*v_j1 - u_j1*v_k2 - u_j2*v_k1 + u_k1*v_j2 + u_k2*v_j1 + u_j2*v_k2 - u_k2*v_j2)^2, -((u_j1 - u_j2)*(u_k1 - u_k2)*(u_j2*v_k1 - u_k1*v_j2 - u_j2*v_k2 + u_k2*v_j2 + u_k1*v_k2 - u_k2*v_k1))/(u_j1*v_k1 - u_k1*v_j1 - u_j1*v_k2 - u_j2*v_k1 + u_k1*v_j2 + u_k2*v_j1 + u_j2*v_k2 - u_k2*v_j2)^2, -((u_k1 - u_k2)*(v_j1 - v_j2)*(u_j1*v_k1 - u_k1*v_j1 - u_j1*v_k2 + u_k2*v_j1 + u_k1*v_k2 - u_k2*v_k1))/(u_j1*v_k1 - u_k1*v_j1 - u_j1*v_k2 - u_j2*v_k1 + u_k1*v_j2 + u_k2*v_j1 + u_j2*v_k2 - u_k2*v_j2)^2, ((u_j1 - u_j2)*(u_k1 - u_k2)*(u_j1*v_k1 - u_k1*v_j1 - u_j1*v_k2 + u_k2*v_j1 + u_k1*v_k2 - u_k2*v_k1))/(u_j1*v_k1 - u_k1*v_j1 - u_j1*v_k2 - u_j2*v_k1 + u_k1*v_j2 + u_k2*v_j1 + u_j2*v_k2 - u_k2*v_j2)^2, ((u_j1 - u_j2)*(v_k1 - v_k2)*(u_j1*v_j2 - u_j2*v_j1 - u_j1*v_k2 + u_k2*v_j1 + u_j2*v_k2 - u_k2*v_j2))/(u_j1*v_k1 - u_k1*v_j1 - u_j1*v_k2 - u_j2*v_k1 + u_k1*v_j2 + u_k2*v_j1 + u_j2*v_k2 - u_k2*v_j2)^2, -((u_j1 - u_j2)*(u_k1 - u_k2)*(u_j1*v_j2 - u_j2*v_j1 - u_j1*v_k2 + u_k2*v_j1 + u_j2*v_k2 - u_k2*v_j2))/(u_j1*v_k1 - u_k1*v_j1 - u_j1*v_k2 - u_j2*v_k1 + u_k1*v_j2 + u_k2*v_j1 + u_j2*v_k2 - u_k2*v_j2)^2, -((u_j1 - u_j2)*(v_k1 - v_k2)*(u_j1*v_j2 - u_j2*v_j1 - u_j1*v_k1 + u_k1*v_j1 + u_j2*v_k1 - u_k1*v_j2))/(u_j1*v_k1 - u_k1*v_j1 - u_j1*v_k2 - u_j2*v_k1 + u_k1*v_j2 + u_k2*v_j1 + u_j2*v_k2 - u_k2*v_j2)^2, ((u_j1 - u_j2)*(u_k1 - u_k2)*(u_j1*v_j2 - u_j2*v_j1 - u_j1*v_k1 + u_k1*v_j1 + u_j2*v_k1 - u_k1*v_j2))/(u_j1*v_k1 - u_k1*v_j1 - u_j1*v_k2 - u_j2*v_k1 + u_k1*v_j2 + u_k2*v_j1 + u_j2*v_k2 - u_k2*v_j2)^2;
((v_j1 - v_j2)*(v_k1 - v_k2)*(u_j2*v_k1 - u_k1*v_j2 - u_j2*v_k2 + u_k2*v_j2 + u_k1*v_k2 - u_k2*v_k1))/(u_j1*v_k1 - u_k1*v_j1 - u_j1*v_k2 - u_j2*v_k1 + u_k1*v_j2 + u_k2*v_j1 + u_j2*v_k2 - u_k2*v_j2)^2, -((u_j1 - u_j2)*(v_k1 - v_k2)*(u_j2*v_k1 - u_k1*v_j2 - u_j2*v_k2 + u_k2*v_j2 + u_k1*v_k2 - u_k2*v_k1))/(u_j1*v_k1 - u_k1*v_j1 - u_j1*v_k2 - u_j2*v_k1 + u_k1*v_j2 + u_k2*v_j1 + u_j2*v_k2 - u_k2*v_j2)^2, -((v_j1 - v_j2)*(v_k1 - v_k2)*(u_j1*v_k1 - u_k1*v_j1 - u_j1*v_k2 + u_k2*v_j1 + u_k1*v_k2 - u_k2*v_k1))/(u_j1*v_k1 - u_k1*v_j1 - u_j1*v_k2 - u_j2*v_k1 + u_k1*v_j2 + u_k2*v_j1 + u_j2*v_k2 - u_k2*v_j2)^2, ((u_j1 - u_j2)*(v_k1 - v_k2)*(u_j1*v_k1 - u_k1*v_j1 - u_j1*v_k2 + u_k2*v_j1 + u_k1*v_k2 - u_k2*v_k1))/(u_j1*v_k1 - u_k1*v_j1 - u_j1*v_k2 - u_j2*v_k1 + u_k1*v_j2 + u_k2*v_j1 + u_j2*v_k2 - u_k2*v_j2)^2, ((v_j1 - v_j2)*(v_k1 - v_k2)*(u_j1*v_j2 - u_j2*v_j1 - u_j1*v_k2 + u_k2*v_j1 + u_j2*v_k2 - u_k2*v_j2))/(u_j1*v_k1 - u_k1*v_j1 - u_j1*v_k2 - u_j2*v_k1 + u_k1*v_j2 + u_k2*v_j1 + u_j2*v_k2 - u_k2*v_j2)^2, -((u_k1 - u_k2)*(v_j1 - v_j2)*(u_j1*v_j2 - u_j2*v_j1 - u_j1*v_k2 + u_k2*v_j1 + u_j2*v_k2 - u_k2*v_j2))/(u_j1*v_k1 - u_k1*v_j1 - u_j1*v_k2 - u_j2*v_k1 + u_k1*v_j2 + u_k2*v_j1 + u_j2*v_k2 - u_k2*v_j2)^2, -((v_j1 - v_j2)*(v_k1 - v_k2)*(u_j1*v_j2 - u_j2*v_j1 - u_j1*v_k1 + u_k1*v_j1 + u_j2*v_k1 - u_k1*v_j2))/(u_j1*v_k1 - u_k1*v_j1 - u_j1*v_k2 - u_j2*v_k1 + u_k1*v_j2 + u_k2*v_j1 + u_j2*v_k2 - u_k2*v_j2)^2, ((u_k1 - u_k2)*(v_j1 - v_j2)*(u_j1*v_j2 - u_j2*v_j1 - u_j1*v_k1 + u_k1*v_j1 + u_j2*v_k1 - u_k1*v_j2))/(u_j1*v_k1 - u_k1*v_j1 - u_j1*v_k2 - u_j2*v_k1 + u_k1*v_j2 + u_k2*v_j1 + u_j2*v_k2 - u_k2*v_j2)^2 ];


