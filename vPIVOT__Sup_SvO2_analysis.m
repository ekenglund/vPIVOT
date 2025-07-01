function vPIVOT__Sup_SvO2_analysis
%%
global vars

% Create images from k-space data

PreRef_ksp = vars.PC.ksp.PreRef; % dimensions are read, phase, echo, repeat, channel, VENC
PreRef_ksp = mean(PreRef_ksp,6);
PreRef_img = fftshift(fft2(ifftshift(PreRef_ksp)));
PreRef_img = squeeze(mean(PreRef_img,4));
PreRef_mag = squeeze(sqrt(sum(abs(PreRef_img).^2,4)));
figure, imagesc([PreRef_mag(:,:,1),PreRef_mag(:,:,2)]), axis image
PreRef_rawphasediff = -squeeze(angle(sum(PreRef_img(:,:,1,:).*conj(PreRef_img(:,:,2,:)),4)));
figure, imagesc(PreRef_rawphasediff), axis image, caxis([-pi/2 pi/2])

PostRef_ksp = vars.PC.ksp.PostRef; % dimensions are read, phase, echo, repeat, channel, VENC
PostRef_ksp = mean(PostRef_ksp,6);
PostRef_img = fftshift(fft2(ifftshift(PostRef_ksp)));
PostRef_img = squeeze(mean(PostRef_img,4));
PostRef_mag = squeeze(sqrt(sum(abs(PostRef_img).^2,4)));
figure, imagesc([PostRef_mag(:,:,1),PostRef_mag(:,:,2)]), axis image
PostRef_rawphasediff = -squeeze(angle(sum(PostRef_img(:,:,1,:).*conj(PostRef_img(:,:,2,:)),4)));
figure, imagesc(PostRef_rawphasediff), axis image, caxis([-pi/2 pi/2])

%%
Data_wPreRef_ksp = vars.PC.ksp.Data_wPreRef; % dimensions are read, phase, echo, repeat, channel, VENC
Data_wPreRef_ksp = mean(Data_wPreRef_ksp,6);
Data_wPreRef_img = fftshift(fft2(ifftshift(Data_wPreRef_ksp)));
Data_wPreRef_mag = squeeze(sqrt(sum(abs(Data_wPreRef_img).^2,5)));
Data_wPreRef_rawphasediff = -squeeze(angle(sum(Data_wPreRef_img(:,:,1,:,:).*conj(Data_wPreRef_img(:,:,2,:,:)),5)));
figure, for k = 1:size(Data_wPreRef_mag,4), subplot(2,1,1), imagesc([Data_wPreRef_mag(:,:,1,k),Data_wPreRef_mag(:,:,2,k)]), axis image, title(k)
subplot(2,1,2), imagesc(Data_wPreRef_rawphasediff(:,:,k)), axis image, caxis([-pi/2 pi/2]), title(k), pause(0.1), drawnow, end

%%

Data_wPostRef_ksp = vars.PC.ksp.Data_wPostRef; % dimensions are read, phase, echo, repeat, channel, VENC
Data_wPostRef_ksp = mean(Data_wPostRef_ksp,6);
Data_wPostRef_img = fftshift(fft2(ifftshift(Data_wPostRef_ksp)));
Data_wPostRef_mag = squeeze(sqrt(sum(abs(Data_wPostRef_img).^2,5)));
figure, for k = 1:size(Data_wPostRef_mag,4), imagesc([Data_wPostRef_mag(:,:,1,k),Data_wPostRef_mag(:,:,2,k)]), axis image, title(k), pause(0.1), drawnow, end
Data_wPostRef_rawphasediff = -squeeze(angle(sum(Data_wPostRef_img(:,:,1,:,:).*conj(Data_wPostRef_img(:,:,2,:,:)),5)));
figure, for k = 1:size(Data_wPostRef_mag,4), imagesc(Data_wPostRef_rawphasediff(:,:,k)), axis image, caxis([-pi/2 pi/2]), title(k), pause(0.1), drawnow, end

%%
Data_mag = Data_wPreRef_mag;
Data_mag(:,:,:,50:end) = Data_wPostRef_mag(:,:,:,50:end);

Data_rawphasediff = Data_wPreRef_rawphasediff;
Data_rawphasediff(:,:,50:end) = Data_wPostRef_rawphasediff(:,:,50:end);

figure, for k = 1:size(Data_mag,4), subplot(2,1,1), imagesc([Data_mag(:,:,1,k),Data_mag(:,:,2,k)]), axis image, title(k)
subplot(2,1,2), imagesc(Data_rawphasediff(:,:,k)), axis image, caxis([-pi/2 pi/2]), title(k), pause(0.1), drawnow, end

