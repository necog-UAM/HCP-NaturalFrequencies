function Omega = find_neighbours(distance)

%% Subs

%% Inside voxels

% N: number of voxels inside the brain. a) 3423 (the 1925 inner voxels have
% all max number of neighbors [19]. b) 3294 (center of the head included)
% c) 1925 (only inner voxels).
% distance: minimal connectivity (1 [7 voxels]),  maximal (1.5 [19 voxels])
    
    cd('Z:\HCP\data\100307\Restin3')
    load source_forward % only for template
    load source_inverse % only for template

voxel_inside = find(source.inside);
    
    for i=1:3
        pos(:,i) = source.pos(:,i)-min(source.pos(:,i))+1;
        dim(i) = max(pos(:,i));
    end
    
    vtempl = zeros(dim);
    dtempl = zeros(size(voxel_inside,1),prod(dim));
    
    for v=1:size(voxel_inside,1)
        temp = zeros(dim);
        p=pos(voxel_inside(v),:);
        vtempl(p(1),p(2),p(3)) = v;
        temp(p(1),p(2),p(3)) = 1;
        dtempl(v,:) = temp(:);
    end


connmat = zeros(size(voxel_inside,1),size(voxel_inside,1));

    for v1=1:size(voxel_inside,1)
        vpos1 = source.pos(voxel_inside(v1),:);
        for v2=1:size(voxel_inside,1)
            vpos2 = source.pos(voxel_inside(v2),:);
            vdif = norm(vpos2-vpos1);
            if vdif<=distance
                connmat(v1,v2) = 1;
            end
        end
    end

connmat = logical(connmat);

    for v=1:size(voxel_inside,1)
        ind = find(vtempl==v);
        [i1,i2,i3] = ind2sub(size(vtempl),ind);
        xx(v)=i1;
        yy(v)=i2;
        zz(v)=i3;
    end

maxconn = max(sum(connmat));

% Inner 1925


% Output
Omega.source_forward = source_forward;
Omega.source = source;
Omega.voxel_inside = voxel_inside;
Omega.dim = dim;
Omega.xx = xx;
Omega.yy = yy;
Omega.zz = zz;
Omega.connmat = connmat;
Omega.dtempl = dtempl;
Omega.maxconn = maxconn;

end

