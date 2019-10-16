%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%

function probs = p_of_lseg_given_v_probability(lsegs,vps,stdev)


%% normalize
ngroups = size(vps,1);
l_length = lsegs(:,7);
lequation = lsegs(:,8:10);
N = size(lsegs,1);

%% probability based consistency measure, by modeling uncertainty of two
%% endpoints of line segment
probs = ones(N,ngroups);
C = 1/(sqrt(2*pi)*stdev);
for j=1:ngroups
    if isstruct(vps(j))
        v = vps(j).vp;
    else
        v = vps(j,:);
    end
    for i =1:N
        L = l_length(i);
        x1 = lsegs(i,1); y1=lsegs(i,3);
        x2 = lsegs(i,2); y2=lsegs(i,4);
        %Y = dist_v_l(lsegs(i,:),v);
        Y = lequation(i,1:3)*v';
        l1 = norm([x1-v(1),y1-v(2)]);
        l2 = norm([x2-v(1),y2-v(2)]); 
        if l1>l2
           X = sqrt(l1^2-Y^2); 
        else
           X = sqrt(l2^2-Y^2); 
        end
        
        if X<=L
            probs(i,j) = 0.00000001;
        else
            probs(i,j) = C*exp( -(L^2*Y^2)/(2*stdev^2 * (X^2+(X-L)^2)) );
        end
    end
    
end


%%
% if isstruct(vps(1))
%     V =  cell2mat({vps.vp}');
%     V = V'
% else
%     V = vps';
% end
% 
% u  =repmat(V(1,:),N,1);
% v  =repmat(V(2,:),N,1);
% u1 =repmat(lsegs(:,1),1,ngroups);
% u2 =repmat(lsegs(:,2),1,ngroups);
% v1 =repmat(lsegs(:,3),1,ngroups);
% v2 =repmat(lsegs(:,4),1,ngroups);
% L = repmat(l_length,1,ngroups);
% the_lambda = (v1-v2).*u + (u2-u1).*v + u1.*v2 - u2.*v1;
% the_gamma = sqrt((u-u1).^2+(v-v1).^2-the_lambda.^2);
% 
% probs = C * exp( -(the_lambda.^2.*L.^2) ./ (2*stdev^2*(2*the_gamma.^2-2*the_gamma.*L.^2+L.^4))  );
% disp('done');


