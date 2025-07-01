function vPIVOT__Sup_PC_analysis
%%
global vars


% Assume first 30 are pre-exercise (use pre for 1-45, then use post)
ksp = zeros(size(vars.PC.ksp.Data_wPreRef)); % dimensions are read, phase, echo, repeat, channel, venc
if size(ksp,4)>45
    ksp(:,:,:,1:45,:,:) = vars.PC.ksp.Data_wPreRef(:,:,:,1:45,:,:);
    ksp(:,:,:,46:end,:,:) = vars.PC.ksp.Data_wPostRef(:,:,:,46:end,:,:);
else
    ksp = vars.PC.ksp.Data_wPreRef;
end

ksp_PreRef = vars.PC.ksp.PreRef;
ksp_PostRef = vars.PC.ksp.PostRef;

ksp_PreRef = squeeze(ksp_PreRef(:,:,1,:,:,:)); % use only the first echo
ksp_PostRef = squeeze(ksp_PostRef(:,:,1,:,:,:)); % use only the first echo

ksp_PreRef = permute(ksp_PreRef,[1,2,5,3,4]); % dimensions are read, phase, venc, repeat, channel to match analysis prior code
ksp_PostRef = permute(ksp_PostRef,[1,2,5,3,4]); % dimensions are read, phase, venc, repeat, channel to match analysis prior code

img_PreRef = fftshift(fft2(ifftshift(ksp_PreRef)));
img_PostRef = fftshift(fft2(ifftshift(ksp_PostRef)));

ksp = squeeze(ksp(:,:,1,:,:,:)); % use only the first echo
ksp = permute(ksp,[1,2,5,3,4]); % dimensions are read, phase, venc, repeat, channel to match analysis prior code
img = fftshift(fft2(ifftshift(ksp)));

venc = 120;
velocity=(venc/pi)*squeeze(angle(sum(img(:,:,1,:,:).*conj(img(:,:,2,:,:)),5)));
magnitude=squeeze(mean(squeeze(sqrt(sum(abs(img).^2,5))),3)); 


velocity_PreRef=(venc/pi)*squeeze(angle(sum(img_PreRef(:,:,1,:,:).*conj(img_PreRef(:,:,2,:,:)),5)));
velocity_PreRef=mean(velocity_PreRef,3);
magnitude_PreRef=squeeze(mean(squeeze(sqrt(sum(abs(img_PreRef).^2,5))),3)); 
magnitude_PreRef=mean(magnitude_PreRef,3);

velocity_PostRef=(venc/pi)*squeeze(angle(sum(img_PostRef(:,:,1,:,:).*conj(img_PostRef(:,:,2,:,:)),5)));
velocity_PostRef=mean(velocity_PostRef,3);
magnitude_PostRef=squeeze(mean(squeeze(sqrt(sum(abs(img_PostRef).^2,5))),3)); 
magnitude_PostRef=mean(magnitude_PostRef,3);

%%
figure, for k = 1:5:size(ksp,4)
    subplot(3,2,1), imagesc(magnitude_PreRef), axis image
    subplot(3,2,2), imagesc(velocity_PreRef), axis image
    subplot(3,2,3), imagesc(magnitude(:,:,k)), axis image, title(k),
    subplot(3,2,4), imagesc(velocity(:,:,k)), axis image, drawnow, pause(0.1)
    subplot(3,2,5), imagesc(magnitude_PostRef), axis image
    subplot(3,2,6), imagesc(velocity_PostRef), axis image
end

%%
figure, imagesc(magnitude_PreRef), axis image
[artery_roi_pre,xa_pre,ya_pre]=roipoly;
[vein_roi2_pre,xv2_pre,yv2_pre]=roipoly;
[vein_roi1_pre,xv1_pre,yv1_pre]=roipoly;


figure, imagesc(magnitude_PostRef), axis image
[artery_roi_post,xa_post,ya_post]=roipoly;
[vein_roi2_post,xv2_post,yv2_post]=roipoly;
[vein_roi1_post,xv1_post,yv1_post]=roipoly;
%%

figure, for k = 1:10:vars.reps
    imagesc(velocity(:,:,k)), axis image, hold on
    if k<45
        plot(xa_pre,ya_pre), plot(xv1_pre,yv1_pre,'k'),plot(xv2_pre,yv2_pre,'r'),hold off, title(k),drawnow, pause(0.1)
    else
        plot(xa_post,ya_post), plot(xv1_post,yv1_post,'k'),plot(xv2_post,yv2_post,'r'),hold off, title(k),drawnow, pause(0.1)
    end

end


% %%
% figure, imagesc(mean(magnitude(:,:,5:25),3)), axis image
% [artery_roi_1,xa_1,ya_1]=roipoly;
% [vein_roi2_1,xv2_1,yv2_1]=roipoly;
% [vein_roi1_1,xv1_1,yv1_1]=roipoly;
% 
% 
% imagesc(mean(magnitude(:,:,75:95),3)), axis image
% hold on, plot(xa_1,ya_1), plot(xv1_1,yv1_1), plot(xv2_1,yv2_1)
% [artery_roi_2,xa_2,ya_2]=roipoly;
% [vein_roi2_2,xv2_2,yv2_2]=roipoly;
% [vein_roi1_2,xv1_2,yv1_2]=roipoly;
% 
% 
% %%
% figure, for k = 1:45
%     imagesc(velocity(:,:,k)), axis image, hold on, caxis([0 120]),
%     plot(xa_1,ya_1), plot(xv1_1,yv1_1,'k'),plot(xv2_1,yv2_1,'r'),hold off, title(k),drawnow, pause(0.1)
% end
% for k = 46:180
%     imagesc(velocity(:,:,k)), axis image, hold on, caxis([0 120]),
%     plot(xa_2,ya_2), plot(xv1_2,yv1_2,'k'),plot(xv2_2,yv2_2,'r'), hold off, title(k),drawnow, pause(0.1)
% end
% 
% 
% for k =1:180
%     if k<46
%         velocity_artery(:,:,k)=velocity(:,:,k).*artery_roi_1;
%         velocity_vein(:,:,k)=velocity(:,:,k).*vein_roi1_1;
%         velocity_vein2(:,:,k)=velocity(:,:,k).*vein_roi2_1;
%     else
%         velocity_artery(:,:,k)=velocity(:,:,k).*artery_roi_2;
%         velocity_vein(:,:,k)=velocity(:,:,k).*vein_roi1_2;
%         velocity_vein2(:,:,k)=velocity(:,:,k).*vein_roi2_2;
%     end
% end

for k =1:45
    velocity_artery(:,:,k)=velocity(:,:,k).*artery_roi_pre;
    velocity_vein(:,:,k)=velocity(:,:,k).*vein_roi1_pre;
    velocity_vein2(:,:,k)=velocity(:,:,k).*vein_roi2_pre;
end
for k =46:vars.reps
    velocity_artery(:,:,k)=velocity(:,:,k).*artery_roi_post;
    velocity_vein(:,:,k)=velocity(:,:,k).*vein_roi1_post;
    velocity_vein2(:,:,k)=velocity(:,:,k).*vein_roi2_post;
end

%
for k =1:vars.reps
    velocity_artery_temp = squeeze(velocity_artery(:,:,k));
    velocity_artery_mean(k) = mean(velocity_artery_temp(find(velocity_artery_temp)));
    velocity_vein_temp = squeeze(velocity_vein(:,:,k));
    velocity_vein_mean(k) = mean(velocity_vein_temp(find(velocity_vein_temp)));
    velocity_vein_temp2 = squeeze(velocity_vein2(:,:,k));
    velocity_vein_mean2(k) = mean(velocity_vein_temp2(find(velocity_vein_temp2)));
   
end

figure, plot(velocity_artery_mean), hold on, plot(velocity_vein_mean),plot(velocity_vein_mean2)


