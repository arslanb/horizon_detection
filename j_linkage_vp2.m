%
% Copyright 2014 by Kitware, Inc. All Rights Reserved. Please refer to
% LICENSE.TXT for licensing information, or contact General Counsel,
% Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.
%
%  Given a set of line segments, find vanishing points using J-linkage

function vps = j_linkage_vp2(im,lsegs,args,n_models)

stdev = args.stdev;
min_consistency = args.min_consistency;

if ~( exist('n_models','var') && ~isempty(n_models) )
   n_models = 3000; 
end
max_rand_samples = 5000;

debug = 0;
[H,W,c]=size(im);

[n_lsegs, nCols] = size(lsegs);

sample_pair_index = randInt(1, n_lsegs, max_rand_samples,  2);
models = zeros(n_models,3);
effective_model_idx = 0;
for i = 1:max_rand_samples
    j = sample_pair_index(i,1);
    k = sample_pair_index(i,2);
    if j == k
        continue;
    end
    cross_vp = cross(lsegs(j,8:10), lsegs(k,8:10));
    cross_vp = cross_vp / cross_vp(3);
    if (cross_vp(1:2)-lsegs(j,[1,3])) * (cross_vp(1:2)-lsegs(j,[2,4]))' > 0 && (cross_vp(1:2)-lsegs(k,[1,3])) * (cross_vp(1:2)-lsegs(k,[2,4]))' > 0 
       effective_model_idx = effective_model_idx + 1;
       models(effective_model_idx,:) = cross_vp;
       if effective_model_idx >=n_models
           break;
       end
    end
end
disp(['Use ', num2str(effective_model_idx), ' effective hypothesized vps']);
models(effective_model_idx+1:end,:) = [];
p = p_of_lseg_given_v_probability(lsegs,models,stdev);


%% start to cluster

% tStart=tic;
M = (p>min_consistency);   
L = diag(true(n_lsegs,1));

members = cell(n_lsegs,2);
for i=1:n_lsegs
   members{i,2} = M(i,:); 
   members{i,1} = i; 
end

intersection_set = zeros(n_lsegs);
union_set = ones(n_lsegs);
for i = 1:n_lsegs
    for j = i+1:n_lsegs
        intersection_set(i,j) = sum( members{i,2} & members{j,2} );
        intersection_set(j,i) = intersection_set(i,j);
        union_set(i,j) = sum( members{i,2} | members{j,2} );
        union_set(j,i) = union_set(i,j);
    end
end

if debug
    hf=figure;
end
% tEnd = toc(tStart)
% tStart=tic;
n_clusters = n_lsegs;
while 1
    d_matrix = (union_set(1:n_clusters,1:n_clusters)-intersection_set(1:n_clusters,1:n_clusters)) ./ union_set(1:n_clusters,1:n_clusters);
    [min_d,ind] = min(d_matrix(:));
    if min_d == 1
        break;
    end
    
    [i,j]=ind2sub(size(d_matrix),ind);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if debug
        imshow(im);
        hold on;
        
        for k = members{i,1}
            plot(lsegs(k,1:2)+W/2,lsegs(k,3:4)+H/2,'r-','LineWidth', 3);
            hold on;
        end
        
        for k = members{j,1}
            plot(lsegs(k,1:2)+W/2,lsegs(k,3:4)+H/2,'g-','LineWidth', 3);
            hold on;
        end
        hold off;
        pause(0.3);
        clf(hf);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    members{i,2} = members{i,2} & members{j,2}; % it is intersect!!
    members{i,1} = [members{i,1},members{j,1}];
    
    
    for k = 1:size(members,1)
        if k~=i && k~=j
            intersection_set(i,k) = sum( members{i,2} & members{k,2} );
            intersection_set(k,i) = intersection_set(i,k);
            union_set(i,k) = sum( members{i,2} | members{k,2} );
            union_set(k,i) = union_set(i,k);
        end
    end
    members(j,:) = members(n_clusters,:);
    intersection_set(j,:) = intersection_set(n_clusters,:);
    intersection_set(:,j) = intersection_set(:,n_clusters);
    
    union_set(j,:) = union_set(n_clusters,:);
    union_set(:,j) = union_set(:,n_clusters);
    
    n_clusters = n_clusters - 1;
end

% tEnd = toc(tStart)
%% construct clusters and vps
for k = 1:n_clusters
   vps(k,1).vp=[];
   vps(k,1).lsegs=lsegs(members{k,1},:);
end





