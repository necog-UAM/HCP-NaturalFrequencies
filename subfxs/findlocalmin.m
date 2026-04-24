function [pks,locs,locmn] = findlocalmin(d,ndim,thr)
%
% Find local minima (only negative peaks) in 1D, 2D or 3D data
% Similar to findlocalmax, but detects valleys instead of peaks
%
% Use as:
%   [pks,locs] = findlocalmin(d,ndim,thr)
%
% Input:
%   d: data
%   ndim: number of dimensions (1:line, 2:plane, 3:volume)
%   thr: threshold for local minima
% 
% Output:
%   pks: local minima
%   locs: indices of local minima
%   locmn: logical matrix indexing local minima
% 
% Almudena Capilla, UAM

if ~exist('thr','var')
    thr = [];
end

if ndim == 1
    if size(d,1) == 1 && size(d,2) > 1
        d = d';
    end
    dx  = diff(squeeze(d(1:end-1,:)),1,1);          
    dx1 = diff(squeeze(d(2:end,:)),1,1);
    sgn = dx.*dx1;
    sgn = [zeros(1,size(sgn,2)); sgn; zeros(1,size(sgn,2))];
    
    d2x = diff(squeeze(d),2,1);
    d2x = [zeros(1,size(d2x,2)); d2x; zeros(1,size(d2x,2))];
    
    if size(d,2) == 1
        pks   = d(sgn(:,1)<0 & d2x(:,1)>0, 1); % only negative peaks (2nd derivative positive)
        locs  = find(sgn(:,1)<0 & d2x(:,1)>0);
        if ~isempty(thr)
            pks2 = pks(pks < thr);
            locs = locs(pks < thr);
            pks = pks2;
        end
        locmn = zeros(size(d));
        locmn(locs) = 1;
        
    else
        locs = {};
        pks  = {};
        locmn = zeros(size(d));
        for i = 1:size(d,2)
            locs{i} = find(sgn(:,i)<0 & d2x(:,i)>0);
            pks{i}  = d(sgn(:,i)<0 & d2x(:,i)>0, i);
            if ~isempty(thr)
                pks2{i} = pks{i}(pks{i} < thr);
                locs{i} = locs{i}(pks{i} < thr);
                pks = pks2;
            end
            if ~isempty(locs{i})
                locmn(locs{i},i) = 1;
            end
        end
        locmn = logical(locmn);
    end
    
elseif ndim == 2
    df  = diff(squeeze(d(1:end-1,:)),1,1);          
    df1 = diff(squeeze(d(2:end,:)),1,1);
    sgn = df.*df1;
    sgn = [zeros(1,size(sgn,2)); sgn; zeros(1,size(sgn,2))];
    
    d2f = diff(squeeze(d),2,1);
    d2f = [zeros(1,size(d2f,2)); d2f; zeros(1,size(d2f,2))];
    
    xlocmn = zeros(size(d));
    for i = 1:size(d,2)
        locs = find(sgn(:,i)<0 & d2f(:,i)>0);
        if ~isempty(locs)
            xlocmn(locs,i) = 1;
        end
    end
    
    df  = diff(squeeze(d(:,1:end-1)),1,2);          
    df1 = diff(squeeze(d(:,2:end)),1,2);
    sgn = df.*df1;
    sgn = [zeros(size(sgn,1),1) sgn zeros(size(sgn,1),1)];
    
    d2f = diff(squeeze(d),2,2);
    d2f = [zeros(size(d2f,1),1) d2f zeros(size(d2f,1),1)];
    
    ylocmn = zeros(size(d));
    for i = 1:size(d,2)
        locs = find(sgn(:,i)<0 & d2f(:,i)>0);
        if ~isempty(locs)
            ylocmn(locs,i) = 1;
        end
    end
    
    locmn = logical(xlocmn.*ylocmn);
    pks   = d(locmn);
    [i,j] = find(locmn==1);
    if ~isempty(thr)
        for p = 1:length(pks)
            if pks(p) > thr
                locmn(i,j) = 0;
            end
        end
        pks   = d(locmn);
        [i,j] = find(locmn==1);
    end
    locs  = [i,j];

elseif ndim == 3
    locmn = zeros(size(d));
    
    for z = 2:size(d,3)-1
        dz = d(:,:,z);
        
        df  = diff(squeeze(dz(1:end-1,:)),1,1);
        df1 = diff(squeeze(dz(2:end,:)),1,1);
        sgn = df.*df1;
        sgn = [zeros(1,size(sgn,2)); sgn; zeros(1,size(sgn,2))];
        
        d2f = diff(squeeze(dz),2,1);
        d2f = [zeros(1,size(d2f,2)); d2f; zeros(1,size(d2f,2))];
        
        xlocmn = zeros(size(dz));
        for i = 1:size(d,2)
            locs = find(sgn(:,i)<0 & d2f(:,i)>0);
            if ~isempty(locs)
                xlocmn(locs,i) = 1;
            end
        end
        
        df  = diff(squeeze(dz(:,1:end-1)),1,2);
        df1 = diff(squeeze(dz(:,2:end)),1,2);
        sgn = df.*df1;
        sgn = [zeros(size(sgn,1),1) sgn zeros(size(sgn,1),1)];
        
        d2f = diff(squeeze(dz),2,2);
        d2f = [zeros(size(d2f,1),1) d2f zeros(size(d2f,1),1)];
        
        ylocmn = zeros(size(dz));
        for i = 1:size(d,2)
            locs = find(sgn(:,i)<0 & d2f(:,i)>0);
            if ~isempty(locs)
                ylocmn(locs,i) = 1;
            end
        end
        
        zlocmn = logical(xlocmn.*ylocmn);
        pks   = dz(zlocmn);
        [i,j] = find(zlocmn==1);
        locs  = [i,j];
        
        for i = 1:length(pks)
            if d(locs(i,1),locs(i,2),z-1) > pks(i) && d(locs(i,1),locs(i,2),z+1) > pks(i)
                locmn(locs(i,1),locs(i,2),z) = 1;
            end
        end
    end
    
    locmn = logical(locmn);
    pks   = d(locmn);
    locs  = [];
    for z = 2:size(d,3)-1
        [i,j] = find(locmn(:,:,z)==1);
        locsz  = [i,j];
        
        if ~isempty(locsz)
            for i = 1:size(locsz,1)
                locs = [locs ; locsz(i,1),locsz(i,2),z];
            end
        end
    end
end
