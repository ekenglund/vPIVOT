function vPIVOT__Inf_T2star_analysis
%%
global vars

% Create images from k-space data
%%
KHdata_ksp = vars.InfGRE.ksp.KHdata(:,2:23,:,:,:); % dimensions are read, phase, echo, repeat, channel
KHdata_img = fftshift(fft2(ifftshift(KHdata_ksp)));
magnitude_KH = squeeze(sqrt(sum(abs(KHdata_img).^2,5)));

ima_PreRef = fftshift(fft2(ifftshift(vars.InfGRE.ksp.PreRef)));
magnitude_PreRef = squeeze(sqrt(sum(abs(ima_PreRef).^2,5)));
magnitude_PreRef = mean(magnitude_PreRef,4);


ima_PostRef = fftshift(fft2(ifftshift(vars.InfGRE.ksp.PostRef)));
magnitude_PostRef = squeeze(sqrt(sum(abs(ima_PostRef).^2,5)));
magnitude_PostRef = mean(magnitude_PostRef,4);

figure, for k = 1:size(magnitude_KH,4)
    imagesc([magnitude_KH(:,:,1,k),magnitude_KH(:,:,2,k),magnitude_KH(:,:,3,k),magnitude_KH(:,:,4,k),magnitude_KH(:,:,5,k)])
    axis image, title(k), drawnow, end

% %%
% for k = 1:size(magnitude_KH,4)-1
%     for j = 1:size(magnitude_KH,3)
%         magnitude_KH_reg(:,:,j,k) = imregtform(magnitude_KH(:,:,j,k), magnitude_KH(:,:,j,end),"rigid");
%     end
% end

%% Pre exercise
figure, subplot(1,2,1), imagesc(magnitude_PreRef(:,:,2)), axis square
subplot(1,2,2), imagesc(mean(magnitude_KH(:,:,2,2:30),4)), axis square
title('Draw ROI around soleus')
[muscle_ROI_pre,x_pre,y_pre] = roipoly;
% Post exercise
figure, subplot(1,2,1), imagesc(magnitude_PostRef(:,:,2)), axis square
subplot(1,2,2), imagesc(mean(magnitude_KH(:,:,2,70:90),4)), axis square
hold on, plot(x_pre,y_pre,'k')
title('Draw ROI around soleus')
[muscle_ROI_post,x_post,y_post] = roipoly;
%%
repeats = size(magnitude_KH,4);
read = size(magnitude_KH,1);
phase = size(magnitude_KH,2);

Signal = zeros(repeats,5);
Signal_e1 = zeros(repeats,1);
Signal_e2 = zeros(repeats,1);
Signal_e3 = zeros(repeats,1);
Signal_e4 = zeros(repeats,1);
Signal_e5 = zeros(repeats,1);

magnitude_KH_e1 = squeeze(magnitude_KH(:,:,1,:));
magnitude_KH_e2 = squeeze(magnitude_KH(:,:,2,:));
magnitude_KH_e3 = squeeze(magnitude_KH(:,:,3,:));
magnitude_KH_e4 = squeeze(magnitude_KH(:,:,4,:));
magnitude_KH_e5 = squeeze(magnitude_KH(:,:,5,:));
for k = 1:45
    Signal_e1(k) = mean(magnitude_KH_e1(find(muscle_ROI_pre)+read*phase*(k-1)));
    Signal_e2(k) = mean(magnitude_KH_e2(find(muscle_ROI_pre)+read*phase*(k-1)));
    Signal_e3(k) = mean(magnitude_KH_e3(find(muscle_ROI_pre)+read*phase*(k-1)));
    Signal_e4(k) = mean(magnitude_KH_e4(find(muscle_ROI_pre)+read*phase*(k-1)));
    Signal_e5(k) = mean(magnitude_KH_e5(find(muscle_ROI_pre)+read*phase*(k-1)));
end
for k = 46:repeats
    Signal_e1(k) = mean(magnitude_KH_e1(find(muscle_ROI_post)+read*phase*(k-1)));
    Signal_e2(k) = mean(magnitude_KH_e2(find(muscle_ROI_post)+read*phase*(k-1)));
    Signal_e3(k) = mean(magnitude_KH_e3(find(muscle_ROI_post)+read*phase*(k-1)));
    Signal_e4(k) = mean(magnitude_KH_e4(find(muscle_ROI_post)+read*phase*(k-1)));
    Signal_e5(k) = mean(magnitude_KH_e5(find(muscle_ROI_post)+read*phase*(k-1)));
end
%%
figure, for k = 1:5:repeats
    imagesc(squeeze(magnitude_KH(:,:,2,k))), caxis([min(min(min(magnitude_KH(:,:,:,2)))), max(max(max(magnitude_KH(:,:,:,2))))]), hold on,
    if k<45
        plot(x_pre,y_pre), hold off, title(k), pause(0.1), drawnow
    else
        plot(x_post,y_post), hold off, title(k), pause(0.1), drawnow
    end
end
clear magnitude_KH_e1 magnitude_KH_e2 magnitude_KH_e3 magnitude_KH_e4 magnitude_KH_e5
Signal(:,1) = Signal_e1;
Signal(:,2) = Signal_e2;
Signal(:,3) = Signal_e3;
Signal(:,4) = Signal_e4;
Signal(:,5) = Signal_e5;
% figure, plot(Signal)
clear Signal_e1 Signal_e2 Signal_e3 Signal_e4 Signal_e5

%%

% figure, plot(TEs,Signal(repeats/2,:),'.')
SignalLN = log(Signal);
% TEs = [3.785, 6.995, 12.325, 19.325, 26.325];
TEs = [3.8, 7.1, 10.4, 17.4, 24.4];

SignalLN1 = SignalLN(:,2:end);
TEs1 = TEs(2:end);

Fit = zeros(repeats,2);
% Error = zeros(repeats);
clear s
for k = 1:repeats
    Fit(k,:) = polyfit(TEs,SignalLN(k,:),1);
    s(k) = regstats(SignalLN(k,:),TEs,'linear');
    yhat(:,k) = s(k).yhat;
    rsq(k) = s(k).rsquare;
end
clear s1
for k = 1:repeats
    Fit1(k,:) = polyfit(TEs1,SignalLN1(k,:),1);
    s1(k) = regstats(SignalLN1(k,:),TEs1,'linear');
    yhat1(:,k) = s1(k).yhat;
    rsq1(k) = s1(k).rsquare;
end

clear T2star
% figure, plot(-Fit(:,1),'.')
T2star = -1./Fit(:,1);
T2star1 = -1./Fit1(:,1);
% figure, plot(2:2:repeats*2,T2star,'.')
% hold on, plot(2:2:repeats*2,T2star1,'r.')
if repeats >20
figure, plot(2:2:repeats*2,T2star/mean(T2star(2:30)),'.'), ginput(1);
hold on, plot(2:2:repeats*2,T2star1/mean(T2star1(2:30)),'r.')
end
%%
rsq_mean = mean(rsq);
if repeats>20
T2star_baseline = mean(T2star(2:30));
end
T2star1Norm = T2star1./mean(T2star1(2:30))*100;
T2starNorm = T2star./mean(T2star(2:30))*100;
%%
vars.InfGRE.T2star.mask_pre = muscle_ROI_pre;
vars.InfGRE.T2star.mask_post = muscle_ROI_post;
vars.InfGRE.T2star.img = magnitude_KH;
vars.InfGRE.T2star.signal = Signal;
vars.InfGRE.T2star.TEs = TEs;
vars.InfGRE.T2star.AverageTimecourse.TEs1_5 = T2star;
vars.InfGRE.T2star.AverageTimecourse.TEs2_5 = T2star1;
vars.InfGRE.T2star.AverageTimecourse.NormTEs1_5 = T2starNorm;
vars.InfGRE.T2star.AverageTimecourse.NormTEs2_5 = T2star1Norm;