%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%

function [hf] = draw_lines_in_img(im,h)

[H,W,c]=size(im);

color = ['r','g','b','m','y','c'];

hf = figure;
imshow(im); hold on;

%% computed horizion
a = h(1);
b = h(2);
c = h(3);

hx1=-W/2; hy1=(-a*hx1-c)/b;
hx2=W/2;  hy2=(-a*hx2-c)/b;

plot([hx1,hx2]+W/2, [hy1,hy2]+H/2,'Color', [0.5, 0, 0.5], 'LineWidth', 7); hold on;
hold off;


