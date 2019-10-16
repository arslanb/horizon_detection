%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%
% Compute horizon given vanishing points and image information
% Input:    vps: data struture for all detected vanishing points
%           imsize: [im_H, im_W]
%           pp: principle point
%           est_f_range:  estiamted focal length range
%           opt_option: choice of algorithm variations: 'min-error' or 'max-log'
% Output:   h:  horizon line in homogeneous cooridante system.  Empty when
% no horizon is detected
%           zenith: the vanishing point(s) that should be the zenith
%           v_vps: data structure for all vertical vanishing points
%           h_vps: data structure for all horizontal vanishing points 
% Authr: Yiliang Xu, Kitware

function [h,zenith,v_vps,h_vps] = horizon(vps,imsize,pp,est_f_range,opt_option)




h=[]; 
H=imsize(1);
W=imsize(2);
picSize = max(W,H);

if isempty(pp)
   pp = [0,0]; 
end
if ~(~isempty(est_f_range) && exist('est_f_range','var'))
   est_f_range = picSize*[0.28, 3.8]; 
end

[zenith, v_vps,h_vps_finite,h_vps_infinite] = split_vps_into_VH(vps,imsize,pp,est_f_range);
h_vps = [h_vps_finite; h_vps_infinite];

disp(['Zenith: ', num2str(uint8(~isempty(zenith)))]);
disp(['h_vps_finite: ', num2str(length(h_vps_finite))] );
disp(['h_vps_infinite: ', num2str(length(h_vps_infinite))] );

if isempty(zenith)  % no zenith
    if length(h_vps) > 1 % there are 2+ horizonal vp, weighted least square
        disp('No zenith, use all 2+ h_vp by weighted least square');
        w=zeros(length(h_vps));
        hV=zeros(length(h_vps),3);
        for i=1:length(h_vps)
            if strcmp(char(opt_option), 'max-log')
                %w(i,i)=1/h_vps(i).mean_gauss_error;
                w(i,i)=size(h_vps(i).lsegs,1);
            elseif strcmp(char(opt_option), 'min-error')
                w(i,i)=1/h_vps(i).tr_pt;
            else
                w(i,i)=1/h_vps(i).tr_pt;
            end
            hV(i,:) = h_vps(i).vp;
        end
        [eigH, lambda] = eig(hV'*w'*w*hV);
        [tmp, smallest] = min(diag(lambda));
        h = eigH(:, smallest);
        h = h / norm(h(1:2));
    else % have 1- horizontal vp
        disp('No zenith, less than 2 h_vp. Return empty horizon');
        h = [];
    end
else % there is zenith
    ver_line = cross(zenith,[pp,1]);
    ver_line = ver_line / norm(ver_line(1:2));
    if ~isempty(h_vps_infinite) % if there is an infinite h_vp
        if length(h_vps_finite) < 1 % there is no finite h_vp
            disp('Has zenith, no finite h_vp, use all infinite_h_vp by weighted least square');
            w=zeros(length(h_vps));
            hV=zeros(length(h_vps),3);
            for i=1:length(h_vps)
                if strcmp(char(opt_option), 'max-log')
                    %w(i,i)=1/h_vps(i).mean_gauss_error;
                    w(i,i)=size(h_vps(i).lsegs,1);
                elseif strcmp(char(opt_option), 'min-error')
                    w(i,i)=1/h_vps(i).tr_pt;
                else
                    w(i,i)=1/h_vps(i).tr_pt;
                end
                hV(i,:) = h_vps(i).vp;
            end
            [eigH, lambda] = eig(hV'*w'*w*hV);
            [tmp, smallest] = min(diag(lambda));
            h = eigH(:, smallest);
            h = h / norm(h(1:2));
        else % there is finite (and zenith), then ignore infinite_h_vp
            % given ver_line is [a,b,c], then horizon should be [-b,a,cc]
            disp('Has zenith, has finite h_vp, ignore infinite h_vp');
            weight_sum = 0;
            cc = 0;
            for j=1:length(h_vps_finite)
                if strcmp(char(opt_option), 'max-log')
                    %weight = 1 / h_vps_finite(j).mean_gauss_error;
                    weight = size(h_vps_finite(j).lsegs,1);
                elseif strcmp(char(opt_option), 'min-error')
                    weight = 1 / h_vps_finite(j).tr_pt;
                else
                    weight = 1 / h_h_vps_finitevps(j).tr_pt;
                end
                weight_sum = weight_sum + weight;
                cc = cc + weight*( -ver_line(1)*h_vps_finite(j).vp(2) + ver_line(2)*h_vps_finite(j).vp(1) );
            end
            cc = cc/weight_sum;
            h = [-ver_line(2), ver_line(1), cc];
            h = h / norm(h(1:2));
        end
    else % there is only finite horizontal vp, there is zenith
        if length(h_vps)>0
            disp('Has zenith, has finite h_vp, no infinite h_vp');
            % given ver_line is [a,b,c], then horizon should be [-b,a,cc]
            weight_sum = 0;
            cc = 0;
            for j=1:length(h_vps)
                if strcmp(char(opt_option), 'max-log')
                    %weight = 1 / h_vps(j).mean_gauss_error;
                    weight = size(h_vps(j).lsegs,1);
                elseif strcmp(char(opt_option), 'min-error')
                    weight = 1 / h_vps(j).tr_pt;
                else
                    weight = 1 / h_vps(j).tr_pt;
                end
                weight_sum = weight_sum + weight;
                cc = cc + weight*( -ver_line(1)*h_vps(j).vp(2) + ver_line(2)*h_vps(j).vp(1) );
            end
            cc = cc/weight_sum;
            h = [-ver_line(2), ver_line(1), cc];
            h = h / norm(h(1:2));
        else
            disp('Has zenith, but no h_vp. Return empty horizon');
            h = [];
        end
    end
end
