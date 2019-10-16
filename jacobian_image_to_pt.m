%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%

function J = jacobian_image_to_pt(varargin)

if nargin == 2
    u = varargin{1}(1);
    v = varargin{1}(2);
    f = varargin{2};
elseif nargin == 3
    u = varargin{1};
    v = varargin{2};
    f = varargin{3};    
else
    error('Incorrect number of inputs for jacobian_image_to_pt');
end

temp1 = f^2+u^2;
temp2 = temp1+v^2;

J = [f/temp1,                   0;
     u*v/(temp2*sqrt(temp1)),  -sqrt(temp1)/temp2];