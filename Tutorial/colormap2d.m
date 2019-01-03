function map2d=colormap2d(nrb);

rbmap = zeros(nrb*8-1,3);
rbmap(1:nrb,3) = (0.5+1/2/nrb:1/nrb/2:1)';
rbmap(nrb+1:4*nrb+1,3) = 1;
rbmap(nrb+1:3*nrb,2) = (1/2/nrb:1/nrb/2:1)';
rbmap(3*nrb+1:4*nrb,2) = 1;
rbmap(3*nrb+1:4*nrb,1) = rbmap(1:nrb,3);
rbmap(4*nrb+1:8*nrb-1,:)=fliplr(flipud(rbmap(1:4*nrb-1,:)));

map2d=rbmap;
