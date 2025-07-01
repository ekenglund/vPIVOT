function backgroundphase=compute_field(rawphasediff,maskedmagnitude)
phi_vec=rawphasediff(:);
magimg_vec=maskedmagnitude(:);

[gX,gY]=ndgrid(1:size(rawphasediff,1),1:size(rawphasediff,2));
gXvec=gX(:);
gYvec=gY(:);

N=length(phi_vec);

fitterms=zeros(N,6);

fitterms(:,1)=ones(N,1);
fitterms(:,2)=gXvec;
fitterms(:,3)=gYvec;
fitterms(:,4)=gXvec.*gYvec;
fitterms(:,5)=gXvec.^2;
fitterms(:,6)=gYvec.^2;

Pfit=find_best_fit(phi_vec,magimg_vec,fitterms);
backgroundphase=zeros(size(rawphasediff));
backgroundphase(:)=Pfit;