%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%

function [zenith,v_vps,h_vps_finite,h_vps_infinite]=split_vps_into_VH(vps,imsize,pp,est_f_range)
v_vps_idx=[];
h_vps_idx_finite=[];
h_vps_idx_infinite=[];
H=imsize(1);
W=imsize(2);
picSize = max(H,W);
for i =1:length(vps)
    if abs(vps(i).vp(2)) >  2*min(W,H)  && abs(vps(i).vp(1)-pp(1))/abs(vps(i).vp(2)-pp(2)) < tan(20/180*pi)
        v_vps_idx(end+1) = i;
    else
        if abs(vps(i).vp(2)-pp(2))/abs(vps(i).vp(1)-pp(1)) < tan(45/180*pi) || ...
                (  vps(i).vp(1) > -W/2 && vps(i).vp(1) < W/2 && vps(i).vp(2) > -0.5*H && vps(i).vp(2) < 0.5*H )
            if abs( vps(i).vp(1) ) > 10*picSize
                h_vps_idx_infinite(end+1) = i;
            else
                h_vps_idx_finite(end+1) = i;
            end
        end
    end
end


v_vps = vps(v_vps_idx);
h_vps_finite = vps(h_vps_idx_finite);
h_vps_infinite = vps(h_vps_idx_infinite);


%% zenith: just pick the largest one
if length(v_vps)<1
    zenith=[];
    v_vps =[];
    h_vps_finite = [h_vps_finite;h_vps_infinite];
else
    largest_idx = 0;
    largest_size = -1;
    for i=1:length(v_vps)
        if largest_size < size( v_vps(i).lsegs, 1 )
            largest_idx = i;
            largest_size = size( v_vps(i).lsegs, 1 );
        end
    end
    zenith = v_vps(largest_idx).vp;
    v_vps=[v_vps(largest_idx)];
    
    
    % filter out finite h_vps based on focal length range
    remove=[];
    ver_line = cross(zenith,[pp,1]);
    ver_line = ver_line / norm(ver_line(1:2));
    for i = 1:length(h_vps_finite)
        c = ( -ver_line(1)*h_vps_finite(i).vp(2) + ver_line(2)*h_vps_finite(i).vp(1) );
        h = [-ver_line(2), ver_line(1), c];
        h = h / norm(h(1:2));
        intersec = cross(ver_line,h);
        intersec = intersec / intersec(3);
        if norm(intersec(1:2)-pp) > 0.07*picSize % consider only when the intersec is not very close to pp
            if (zenith(1:2)-pp)*(intersec(1:2)-pp)' > 0
                remove(end+1) = i;
            else
                f = sqrt( norm(zenith(1:2)-pp) * norm(intersec(1:2)-pp) );
                if f<est_f_range(1) || f>est_f_range(2)
                    remove(end+1)=i;
                end
            end
        end
    end
    
    h_vps_finite(remove) = [];
end



