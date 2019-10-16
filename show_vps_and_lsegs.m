%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%

function fh=show_vps_and_lsegs(vps,im,color)

if isempty(color)
   color = [0,0,0;
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
end



fh=figure;
n_vps = length(vps);
[H,W,c]=size(im);
colorsize = size(color,1);


% for i=1:n_vps
% %     if ~isempty(vps(i).vp)
% %         plot(vps(i).vp(1)*H+W/2,vps(i).vp(2)*H+H/2,[color(mod(i,colorsize)+1),'*']);
% %         hold on;
% %     end
%     n_lsegs = size(vps(i).lsegs,1);
%     for j = 1:n_lsegs
%         plot(vps(i).lsegs(j,1:2)*H+W/2,vps(i).lsegs(j,3:4)*H+H/2,[color(mod(i,colorsize)+1),'-'],'LineWidth', 1);
%         hold on;
%     end
% end


imshow(im);
hold on;

for i=1:n_vps
%     if ~isempty(vps(i).vp)
%         plot(vps(i).vp(1)*H+W/2,vps(i).vp(2)*H+H/2,[color(mod(i,colorsize)+1),'*']);
%         hold on;
%     end
    n_lsegs = size(vps(i).lsegs,1);
    for j = 1:n_lsegs
        plot(vps(i).lsegs(j,1:2)+W/2,vps(i).lsegs(j,3:4)+H/2,'Color',color(mod(i,colorsize)+1,:),'LineWidth', 2);
        hold on;
    end
end

hold off;






