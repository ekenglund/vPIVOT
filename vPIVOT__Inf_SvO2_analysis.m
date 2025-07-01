function vPIVOT__Inf_SvO2_analysis
%%
global vars

% Create images from k-space data
%%
% Assume first 30 are pre-exercise (use pre for 1-45, then use post)
ksp = zeros(size(vars.InfGRE.ksp.Data_wPreRef)); % dimensions are read, phase, echo, repeat, channel
repeats = size(ksp,4);
ksp(:,:,:,1:45,:) = vars.InfGRE.ksp.Data_wPreRef(:,:,:,1:45,:);
ksp(:,:,:,46:end,:) = vars.InfGRE.ksp.Data_wPostRef(:,:,:,46:end,:);
% ksp = vars.InfGRE.ksp.Data_wPostRef;

ima_PreRef = fftshift(fft2(ifftshift(vars.InfGRE.ksp.PreRef)));
ima_PostRef = fftshift(fft2(ifftshift(vars.InfGRE.ksp.PostRef)));
ima = fftshift(fft2(ifftshift(ksp)));

magnitude_dyn = squeeze(sqrt(sum(abs(ima).^2,5)));
magnitude_PreRef = squeeze(sqrt(sum(abs(ima_PreRef).^2,5)));
magnitude_PreRef2 = mean(magnitude_PreRef,4);
magnitude_PostRef = squeeze(sqrt(sum(abs(ima_PostRef).^2,5))); % sum of squares coil combo
magnitude_PostRef2 = mean(magnitude_PostRef,4);

magnitude(:,:,:,1:3) = magnitude_PreRef;
magnitude(:,:,:,4:repeats+3) = magnitude_dyn;
magnitude(:,:,:,repeats+4:repeats+6) = magnitude_PostRef;

rawphasediff_dyn = -squeeze(angle(sum(ima(:,:,1,:,:).*conj(ima(:,:,2,:,:)),5)));
rawphasediff_PreRef = -squeeze(angle(sum(ima_PreRef(:,:,1,:,:).*conj(ima_PreRef(:,:,2,:,:)),5)));
rawphasediff_PostRef = -squeeze(angle(sum(ima_PostRef(:,:,1,:,:).*conj(ima_PostRef(:,:,2,:,:)),5)));

rawphasediff(:,:,1:3) = rawphasediff_PreRef;
rawphasediff(:,:,4:repeats+3) = rawphasediff_dyn;
rawphasediff(:,:,repeats+4:repeats+6) = rawphasediff_PostRef;

%%
figure, imagesc(magnitude(:,:,1,100)), axis square, colormap gray,
disp(' create background tissue mask');
figure, imagesc([rawphasediff(:,:,10), rawphasediff(:,:,60), rawphasediff(:,:,100)]), axis image

unwrap = input('need to unwrap?  (1=yes)');


if unwrap == 1
    title('select region that is wrapped')
    unwrap_roi = roipoly;

    rawphasediff_unwrap = rawphasediff; 
    rawphasediff_wrap = rawphasediff;

    unwrap_roi = unwrap_roi(1:96,1:96);

%     rawphasediff_PreRef_unwrap = rawphasediff_PreRef;
%     rawphasediff_PreRef_wrap = rawphasediff_PreRef;
% 
%     rawphasediff_PostRef_unwrap = rawphasediff_PostRef;
%     rawphasediff_PostRef_wrap = rawphasediff_PostRef;
    for j = 1:186
        for kr = 1:size(unwrap_roi,1)
            for kc = 1:size(unwrap_roi,2)
                if rawphasediff_wrap(kr,kc,j)<-1
                    rawphasediff_unwrap(kr,kc,j) = rawphasediff_wrap(kr,kc,j)+2*pi;
%                 if j == 1
%                     if rawphasediff_PreRef(kr,kc,j)<-1
%                     rawphasediff_PreRef_unwrap(kr,kc,j) = rawphasediff_PreRef(kr,kc,j)+2*pi;
%                     end
%                     if rawphasediff_PostRef(kr,kc,j)<-1
%                     rawphasediff_PostRef_unwrap(kr,kc,j) = rawphasediff_PostRef(kr,kc,j)+2*pi;
%                     end
%                 end
                end
            end
        end
        end
    rawphasediff = rawphasediff_unwrap;
%     rawphasediff_PreRef = rawphasediff_PreRef_unwrap;
%     rawphasediff_PostRef = rawphasediff_PostRef_unwrap;
end
% rawphasediff now contains the reference images (first three and last
% three)
%%
basemask=mean(rawphasediff(:,:,repeats+4:repeats+6),3);
figure, imagesc(basemask), axis square, caxis([min(rawphasediff(:))+1, max(rawphasediff(:))-1])
backgroundouterbinary=roipoly;
%
nVes = input('How many vessels do you want to mask out?  ');
basemask2 = basemask;
close, figure, imagesc(basemask2), caxis([min(rawphasediff(:))+1, max(rawphasediff(:))-1]), axis square
for k = 1:nVes
    title(['Mask out ',num2str(nVes+1-k),' vessels'])
    [y,x]=ginput(1);
    for dx = -(x-1):96-x;
        for dy = -(y-1):96-y;
            dist = sqrt(dx^2+dy^2);
            if dist <=3
                backgroundouterbinary(x+dx,y+dy)=0;
                basemask2(x+dx,y+dy) = 0;          
            end
        end
    end
    imagesc(basemask2), axis image
    hold on, plot(y,x,'.','MarkerSize',25)
    hold off
end
%%
backgroundmask=basemask2.*backgroundouterbinary;
read = 96; phase = 96; repeats = repeats+6;
imshow(backgroundmask);
% write_mda('backgroundmask.mda',backgroundmask);
% clear basemask backgroundmask;
maskedmagnitude=zeros(read,phase,repeats);
for k=1:repeats;
    maskedmagnitude(:,:,k)=backgroundmask.*magnitude(:,:,k);
end;
% clear backgroundbinary k;
% figure, for k = 1:10:repeats
% imagesc(maskedmagnitude(:,:,k)),title(k), drawnow, pause(0.1)
% end

% Field fitting 
disp(' field fitting');
backgroundphase=zeros(read,phase,repeats);
for k=1:repeats
    backgroundphase(:,:,k)=compute_field(rawphasediff(:,:,k),maskedmagnitude(:,:,k));
end;
phasediff=rawphasediff-backgroundphase;
% write_mda('phasediff.mda',phasediff);
figure, for k = 1:10:repeats
    imagesc(phasediff(:,:,k)), axis image, caxis([-pi/2 pi/2]), title(k), drawnow, pause(0.1)
end
% clear maskedmagnitude rawphasediff backgroundphase k;

%% Vessel ROI determination
% disp(' vessel ROI determination');
% close all,clear all,clc
% 
% load('SvO2_Location.mat');
% phasediff = readmda('phasediff.mda');
% magnitude = readmda('magnitude.mda');
% repeats = 186; read = 96; phase = 96;
gamma=42.576*2*pi*1e6;
chioxy=-0.008*4*pi*1e-6;
chido=0.273*4*pi*1e-6;
b0=2.89;
deltaTE=0.00327;
% SvO2 ROI determination
disp('  SvO2 ROI determiantion');
% 
% figure, imagesc(Stack(:,:,2)), 
% h = imellipse(gca,pos(2,:));
figure, subplot(1,2,1),  imagesc(squeeze(mean(magnitude(:,:,1,1:3),4))), axis square
subplot(1,2,2),imagesc(mean(phasediff(:,:,1:3),3)), caxis([-pi/2, pi/2]), axis image, 
% if exist('SvO2_ROIs.mat')==2
%     load('SvO2_ROIs.mat');
% else
title('select vein pre-exercise')
% hold on, plot(xv,yv)
[Vein_pre,xv_pre,yv_pre]=roipoly;
title('select background tissue')
[BackgroundTissue1_pre,xt1_pre,yt1_pre] = roipoly;
[BackgroundTissue2_pre,xt2_pre,yt2_pre] = roipoly;
% end
%%

figure, subplot(1,2,1),  imagesc(squeeze(mean(magnitude(:,:,1,65:70),4))), axis square
subplot(1,2,2),imagesc(mean(phasediff(:,:,65:70),3)), caxis([-pi/2, pi/2]), axis image, 
% if exist('SvO2_ROIs.mat')==2
%     load('SvO2_ROIs.mat');
% else
hold on, plot(xv_pre,yv_pre), plot(xt1_pre,yt1_pre,'k'),plot(xt2_pre,yt2_pre,'k'),
title('select vein post-exercise')
% hold on, plot(xv,yv)
[Vein_post,xv_post,yv_post]=roipoly;
title('select background tissue')
[BackgroundTissue1_post,xt1_post,yt1_post] = roipoly;
[BackgroundTissue2_post,xt2_post,yt2_post] = roipoly;
% end


%%
%
figure, for k = 1:5:45
    imagesc(phasediff(:,:,k)), axis image, hold on
    hold on, plot(xv_pre,yv_pre), plot(xt1_pre,yt1_pre,'k'),plot(xt2_pre,yt2_pre,'k'), caxis([-pi pi]),hold off, title(k),drawnow, pause(0.1)
end
for k = 46:5:repeats
    imagesc(phasediff(:,:,k)), axis image, hold on
    hold on, plot(xv_post,yv_post), plot(xt1_post,yt1_post,'k'),plot(xt2_post,yt2_post,'k'), caxis([-pi pi]),hold off, title(k),drawnow, pause(0.1)
end
%%
for k = 1:45
    BackgroundTissue(:,:,k) = BackgroundTissue1_pre+BackgroundTissue2_pre;
    Vein(:,:,k) = Vein_pre;
end
for k = 46:repeats
    BackgroundTissue(:,:,k) = BackgroundTissue1_post+BackgroundTissue2_post;
    Vein(:,:,k) = Vein_post;
end
%%
for k = 1:repeats
    VeinPhaseIma(:,:,k) = phasediff(:,:,k).*Vein(:,:,k);
    TissuePhaseIma(:,:,k) = phasediff(:,:,k).*BackgroundTissue(:,:,k);
end

for k = 1:repeats
    Vein_temp = Vein(:,:,k);
    BackgroundTissue_temp = BackgroundTissue(:,:,k);
    mean_VeinPhase(k) = mean(VeinPhaseIma(find(Vein_temp)+read*phase*(k-1)));
    mean_TissuePhase(k) = mean(TissuePhaseIma(find(BackgroundTissue_temp)+read*phase*(k-1)));

    median_VeinPhase(k) = median(VeinPhaseIma(find(Vein_temp)+read*phase*(k-1)));
    median_TissuePhase(k) = median(TissuePhaseIma(find(BackgroundTissue_temp)+read*phase*(k-1)));
end

figure, plot(mean_VeinPhase), hold on, plot(mean_TissuePhase,'k')
plot(median_VeinPhase), hold on, plot(median_TissuePhase)
%
deltaPhi = mean_VeinPhase-mean_TissuePhase;
% deltaPhi_smooth = smoothdata(mean_VeinPhase)-smoothdata(mean_TissuePhase);
deltaPhi_median = median_VeinPhase-median_TissuePhase;
% deltaPhi_median_smooth = smoothdata(median_VeinPhase)-smoothdata(median_TissuePhase);
hct = input('Hct? ');
vessel_angle = input('Vessel Angle? ');
SvO2=(1-2*(abs(deltaPhi))/(gamma*chido*b0*deltaTE*((cos(vessel_angle*(pi/180)))^2-1/3)*hct)+chioxy/chido)*100;
SvO2_median=(1-2*(abs(deltaPhi_median))/(gamma*chido*b0*deltaTE*((cos(vessel_angle*(pi/180)))^2-1/3)*hct)+chioxy/chido)*100;
% SvO2_smooth=(1-2*(abs(deltaPhi_smooth))/(gamma*chido*b0*deltaTE*((cos(vessel_angle*(pi/180)))^2-1/3)*hct)+chioxy/chido)*100;
% SvO2_median_smooth=(1-2*(abs(deltaPhi_median_smooth))/(gamma*chido*b0*deltaTE*((cos(vessel_angle*(pi/180)))^2-1/3)*hct)+chioxy/chido)*100;
figure, plot(SvO2), hold on, plot(SvO2_median)%,plot(SvO2_smooth), plot(SvO2_median_smooth)
%%

vars.InfGRE.SvO2.img.magnitude = magnitude;
vars.InfGRE.SvO2.img.rawphasediff = rawphasediff;
vars.InfGRE.SvO2.img.backgroundmask = backgroundmask;
vars.InfGRE.SvO2.img.phasediff = phasediff;
vars.InfGRE.SvO2.mask.BackgroundTissue = BackgroundTissue;
vars.InfGRE.SvO2.mask.Vein = Vein;
vars.InfGRE.SvO2.params.hct = hct;
vars.InfGRE.SvO2.params.vessel_angle = vessel_angle;
vars.InfGRE.SvO2.params.deltaTE = deltaTE;
vars.InfGRE.SvO2.params.gamma = gamma;
vars.InfGRE.SvO2.params.chioxy = chioxy;
vars.InfGRE.SvO2.params.chido = chido;
vars.InfGRE.SvO2.params.b0 = b0;
vars.InfGRE.SvO2.AverageTimecourse.mean_VeinPhase = mean_VeinPhase';
vars.InfGRE.SvO2.AverageTimecourse.mean_TissuePhase = mean_TissuePhase';
vars.InfGRE.SvO2.AverageTimecourse.deltaPhi = deltaPhi';
vars.InfGRE.SvO2.AverageTimecourse.SvO2 = SvO2';
% vars.InfGRE.SvO2.AverageTimecourse.SvO2_smooth = SvO2_smooth';
vars.InfGRE.SvO2.AverageTimecourse.SvO2_median = SvO2_median';
% vars.InfGRE.SvO2.AverageTimecourse.SvO2_median_smooth = SvO2_median_smooth';