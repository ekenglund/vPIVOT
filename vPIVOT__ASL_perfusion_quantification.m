function vPIVOT__ASL_perfusion_quantification
%%
global vars

read = 80; fullphase = 80; channels =size(vars.ASL.raw.EPI,3); repeats = vars.reps;
ImaNS = vars.ASL.img.ImaNS; ImaSS = vars.ASL.img.ImaSS;
PLD = 952.19; lambda = 90; % TE = 7.4;
T = PLD/(1000*60);
T1t = 1420/(1000*60);
deltaM = zeros(read,fullphase,repeats-1);
% motion correction added 1/2/2024
[optimizer, metric] = imregconfig('multimodal');
optimizer.MaximumIterations = 100;
% 
% tested a variety of different strategies - better to first register to
% the original image and then between dynamics? better to register the SS
% and NS together or separate? lots of potential for iterating. At this
% point, for the sake of time, we will register NS and SS separately
% frame-by-frame, then together register to SS final image as fixed. For
% some reason monomodal optimizer didn't work so we will use multimodal.


ImaNS_MC1 = zeros(read,fullphase,repeats); % order is NS then SS
ImaSS_MC1 = zeros(read,fullphase,repeats);
% ImaNS_MC2 = zeros(read,fullphase,repeats); % order is NS then SS
% ImaSS_MC2 = zeros(read,fullphase,repeats);
% ImaNS_MC3 = zeros(read,fullphase,repeats); % order is NS then SS
% ImaSS_MC3 = zeros(read,fullphase,repeats);
% ImaNS_MC4 = zeros(read,fullphase,repeats); % order is NS then SS
% ImaSS_MC4 = zeros(read,fullphase,repeats);
% ImaNS_MC5 = zeros(read,fullphase,repeats); % order is NS then SS
% ImaSS_MC5 = zeros(read,fullphase,repeats);

for k = 1:repeats-1
    ImaNS_MC1(:,:,k) = imregister(ImaNS(:,:,k+1),ImaNS(:,:,k),"rigid",optimizer,metric);
    ImaSS_MC1(:,:,k) = imregister(ImaSS(:,:,k+1),ImaNS(:,:,k),"rigid",optimizer,metric);
end
% 1
% 
% for k = 1:repeats-1
%     ImaNS_MC2(:,:,k) = imregister(ImaNS_MC1(:,:,k+1),ImaNS_MC1(:,:,k),"rigid",optimizer,metric);
%     ImaSS_MC2(:,:,k) = imregister(ImaSS_MC1(:,:,k+1),ImaNS_MC1(:,:,k),"rigid",optimizer,metric);
% end
% 2
% for k = 1:repeats-1
%     ImaNS_MC3(:,:,k) = imregister(ImaNS_MC2(:,:,k+1),ImaNS_MC2(:,:,k),"rigid",optimizer,metric);
%     ImaSS_MC3(:,:,k) = imregister(ImaSS_MC2(:,:,k+1),ImaNS_MC2(:,:,k),"rigid",optimizer,metric);
% end
% 3
% for k = 1:repeats-1
%     ImaNS_MC4(:,:,k) = imregister(ImaNS_MC3(:,:,k+1),ImaNS_MC3(:,:,k),"rigid",optimizer,metric);
%     ImaSS_MC4(:,:,k) = imregister(ImaSS_MC3(:,:,k+1),ImaNS_MC3(:,:,k),"rigid",optimizer,metric);
% end
% 4
% for k = 1:repeats-1
%     ImaNS_MC5(:,:,k) = imregister(ImaNS_MC4(:,:,k+1),ImaNS_MC4(:,:,k),"rigid",optimizer,metric);
%     ImaSS_MC5(:,:,k) = imregister(ImaSS_MC4(:,:,k+1),ImaNS_MC4(:,:,k),"rigid",optimizer,metric);
% end
% 5
% 
% for k = 1:repeats-1
%     ImaNS_MC1_Smoothed(:,:,k) = imgaussfilt(ImaNS_MC1(:,:,k),1);
%     ImaSS_MC1_Smoothed(:,:,k) = imgaussfilt(ImaSS_MC1(:,:,k),1);
%     ImaNS_MC2_Smoothed(:,:,k) = imgaussfilt(ImaNS_MC2(:,:,k),1);
%     ImaSS_MC2_Smoothed(:,:,k) = imgaussfilt(ImaSS_MC2(:,:,k),1);
% end


% 
% 
% 
% ImaNS_MC2 = zeros(read,fullphase,repeats); % order is NS then SS
% ImaSS_MC2 = zeros(read,fullphase,repeats);
% 
% for k = 1:size(ImaNS_MC1,3)
%     ImaNS_MC2(:,:,k) = imregister(ImaNS_MC1(:,:,k),ImaNS_MC1(:,:,1),"rigid",optimizer,metric);
%     ImaSS_MC2(:,:,k) = imregister(ImaSS_MC1(:,:,k),ImaSS_MC1(:,:,1),"rigid",optimizer,metric);
% end


% 
% 
% 
% Ima_fulltimeline(:,:,1:2:repeats*2) = ImaNS;
% Ima_fulltimeline(:,:,2:2:repeats*2) = ImaSS;
% Ima_fulltimeline_fixed = Ima_fulltimeline(:,:,1);
% Ima_fulltimeline_MC = zeros(read,fullphase,repeats*2);
% Ima_fulltimeline_MC(:,:,1) = Ima_fulltimeline(:,:,1);


% for k = 2:repeats
%     ImaNS_MC(:,:,k) = imregister(ImaNS(:,:,k),ImaNS_fixed,"rigid",optimizer,metric);
%     ImaSS_MC(:,:,k) = imregister(ImaSS(:,:,k),ImaSS_fixed,"rigid",optimizer,metric);
% end
% 
% for k = 1:repeats-1
%     ImaNS_MC2(:,:,k) = imregister(ImaNS_MC(:,:,k+1),ImaNS_MC(:,:,k),"rigid",optimizer,metric);
%     ImaSS_MC2(:,:,k) = imregister(ImaSS_MC(:,:,k+1),ImaSS_MC(:,:,k),"rigid",optimizer,metric);
% end
% 
% for k = 2:repeats*2
%     Ima_fulltimeline_MC(:,:,k) = imregister(Ima_fulltimeline(:,:,k),Ima_fulltimeline_fixed,"rigid",optimizer,metric);
% end

% for k = 1:repeats*2-1
%     Ima_fulltimeline_MC2(:,:,k) = imregister(Ima_fulltimeline_MC(:,:,k+1),Ima_fulltimeline_MC(:,:,k),"rigid",optimizer,metric);
% end
%%
ImaSS = ImaSS_MC1;%imgaussfilt(ImaSS_MC1,1);
ImaNS = ImaNS_MC1;%imgaussfilt(ImaNS_MC1,1);
% ImaSS = imgaussfilt(ImaSS,1);
% ImaNS = imgaussfilt(ImaNS,1);
%%
% sumM = zeros(read,fullphase,repeats-1);
% M = zeros(read,fullphase,repeats-1);
clear r p k
for r = 1:read
    for p = 1:fullphase
        for k = 1:repeats-1
            deltaM(r,p,k) = ImaSS(r,p,k)-(ImaNS(r,p,k)+ImaNS(r,p,k+1))/2;
            sumM(r,p,k) = ImaSS(r,p,k)+(ImaNS(r,p,k)+ImaNS(r,p,k+1))/2;
%             deltaM(r,p,k) = ImaSS_MC1(r,p,k)-(ImaNS_MC1(r,p,k)+ImaNS_MC1(r,p,k+1))/2;
%             sumM(r,p,k) = ImaSS_MC1(r,p,k)+(ImaNS_MC1(r,p,k)+ImaNS_MC1(r,p,k+1))/2;
%             deltaM_MC2(r,p,k) = ImaSS_MC2(r,p,k)-(ImaNS_MC2(r,p,k)+ImaNS_MC2(r,p,k+1))/2;
%             sumM_MC2(r,p,k) = ImaSS_MC2(r,p,k)+(ImaNS_MC2(r,p,k)+ImaNS_MC2(r,p,k+1))/2;
        end
    end
end
% load('muscleROIs2.mat')

% figure
% for k = 1:10:repeats-1
%     imagesc(abs(M(:,:,k))), axis image, drawnow, pause(0.1)
% end
% Mmax = max(deltaM(:,:,:end),[],3);
%%
mean_deltaM = mean(deltaM(:,:,60:79),3);
mean_sumM = mean(sumM(:,:,60:79),3);
% figure, imagesc((Mmean<3e-6).*ImaNS(:,:,50)), axis image
% threshold = input('Threshold?  ');
threshold = 100e-6;
deltaM_threshold1 = real(mean_deltaM)<threshold;
deltaM_threshold2 = real(mean_deltaM)>0;
sumM_threshold = real(mean_sumM)>1.3e-5;

Mthreshold = deltaM_threshold1.*deltaM_threshold2.*sumM_threshold;
figure, subplot(2,2,1), imagesc(deltaM_threshold1)
subplot(2,2,2), imagesc(deltaM_threshold2)
subplot(2,2,3), imagesc(sumM_threshold)
subplot(2,2,4), imagesc(Mthreshold)
clear Mthreshold1 Mthreshold2

%%

figure, imagesc(ImaNS(:,:,60).*Mthreshold), axis image, colormap gray

Perf_Img = zeros(size(ImaNS));
for iread = 1:read
    for iphase = 1:fullphase
        for irep = 2:repeats-1
            if sumM_threshold(iread,iphase) == 1;
                Perf_Img(iread,iphase,irep) = -lambda/T*log(((ImaSS(iread,iphase,irep)-(ImaNS(iread,iphase,irep)+ImaNS(iread,iphase,irep+1))/2)/(ImaSS(iread,iphase,irep)+(ImaNS(iread,iphase,irep)+ImaNS(iread,iphase,irep+1))/2))*(1-exp(T/T1t))+1);
                if Perf_Img(iread,iphase,irep)<0
                    Perf_Img(iread,iphase,irep) = NaN;
                elseif Perf_Img(iread,iphase,irep)>250
                    Perf_Img(iread,iphase,irep) = NaN;
                end
            end
        end
    end
end

figure, imagesc(mean(Perf_Img(:,:,60:79),3,'omitnan')), axis image


% if length(dir('*.mat'))==1
%     a = dir('*.mat');
%     load(a.name)
%     mgMask = varsASL.mask.mgMask;
%     lgMask = varsASL.mask.lgMask;
%     gMask = varsASL.mask.gMask;
%     sMask = varsASL.mask.sMask;
%     pMask = varsASL.mask.pMask;
%     tMask = varsASL.mask.tMask;
%     lMask = varsASL.mask.lMask;
%     clear init
%     clear varsASL


% if exist('muscleLEG_ROIs.mat')==2
%     load('muscleLEG_ROIs.mat');
% else
    title('Select ROI around medial gastroc')
    % % hold on, plot(xg,yg)
    [mgMask,xmg,ymg] = roipoly;
    hold on, plot(xmg,ymg)
    title('Select ROI around lateral gastroc')
    % % hold on, plot(xg,yg)
    [lgMask,xlg,ylg] = roipoly;
    hold on, plot(xlg,ylg)
    title('Select ROI around soleus')
    % hold on, plot(xs,ys,'r')
    [sMask,xs,ys] = roipoly;
    hold on, plot(xs,ys,'r')
    title('Select ROI around peroneus')
    % hold on, plot(xp,yp,'r')
    [pMask,xp,yp] = roipoly;
    hold on, plot(xp,yp,'c')
    title('Select ROI around tibialis anterior')
    % hold on, plot(xt,yt,'r')
    [tMask,xt,yt] = roipoly;
    hold on, plot(xt,yt,'g')
% end
gMask = mgMask+lgMask;
lMask = gMask+sMask+pMask+tMask;

%%

mgNSIma = zeros(read,fullphase,repeats);
mgSSIma = zeros(read,fullphase,repeats);
lgNSIma = zeros(read,fullphase,repeats);
lgSSIma = zeros(read,fullphase,repeats);
gNSIma = zeros(read,fullphase,repeats);
gSSIma = zeros(read,fullphase,repeats);
sNSIma = zeros(read,fullphase,repeats);
sSSIma = zeros(read,fullphase,repeats);
pNSIma = zeros(read,fullphase,repeats);
pSSIma = zeros(read,fullphase,repeats);
tNSIma = zeros(read,fullphase,repeats);
tSSIma = zeros(read,fullphase,repeats);
lNSIma = zeros(read,fullphase,repeats);
lSSIma = zeros(read,fullphase,repeats);
NS = zeros(repeats,5);
SS = zeros(repeats,5);
for k = 1:repeats
    mgNSIma(:,:,k) = ImaNS(:,:,k).*mgMask.*Mthreshold;
    NS(k,1) = mean(mgNSIma(find(mgNSIma(:,:,k))+read*fullphase*(k-1)));
    lgNSIma(:,:,k) = ImaNS(:,:,k).*lgMask.*Mthreshold;
    NS(k,2) = mean(lgNSIma(find(lgNSIma(:,:,k))+read*fullphase*(k-1)));
    gNSIma(:,:,k) = ImaNS(:,:,k).*gMask.*Mthreshold;
    NS(k,3) = mean(gNSIma(find(gNSIma(:,:,k))+read*fullphase*(k-1)));
    sNSIma(:,:,k) = ImaNS(:,:,k).*sMask.*Mthreshold;
    NS(k,4) = mean(sNSIma(find(sNSIma(:,:,k))+read*fullphase*(k-1)));
    pNSIma(:,:,k) = ImaNS(:,:,k).*pMask.*Mthreshold;
    NS(k,5) = mean(pNSIma(find(pNSIma(:,:,k))+read*fullphase*(k-1)));
    tNSIma(:,:,k) = ImaNS(:,:,k).*tMask.*Mthreshold;
    NS(k,6) = mean(tNSIma(find(tNSIma(:,:,k))+read*fullphase*(k-1)));
    lNSIma(:,:,k) = ImaNS(:,:,k).*lMask.*Mthreshold;
    NS(k,7) = mean(lNSIma(find(lNSIma(:,:,k))+read*fullphase*(k-1)));
    
    mgSSIma(:,:,k) = ImaSS(:,:,k).*mgMask.*Mthreshold;
    SS(k,1) = mean(mgSSIma(find(mgSSIma(:,:,k))+read*fullphase*(k-1)));
    lgSSIma(:,:,k) = ImaSS(:,:,k).*lgMask.*Mthreshold;
    SS(k,2) = mean(lgSSIma(find(lgSSIma(:,:,k))+read*fullphase*(k-1)));
    gSSIma(:,:,k) = ImaSS(:,:,k).*gMask.*Mthreshold;
    SS(k,3) = mean(gSSIma(find(gSSIma(:,:,k))+read*fullphase*(k-1)));
    sSSIma(:,:,k) = ImaSS(:,:,k).*sMask.*Mthreshold;
    SS(k,4) = mean(sSSIma(find(sSSIma(:,:,k))+read*fullphase*(k-1)));
    pSSIma(:,:,k) = ImaSS(:,:,k).*pMask.*Mthreshold;
    SS(k,5) = mean(pSSIma(find(pSSIma(:,:,k))+read*fullphase*(k-1)));
    tSSIma(:,:,k) = ImaSS(:,:,k).*tMask.*Mthreshold;
    SS(k,6) = mean(tSSIma(find(tSSIma(:,:,k))+read*fullphase*(k-1))); 
    lSSIma(:,:,k) = ImaSS(:,:,k).*lMask.*Mthreshold;
    SS(k,7) = mean(lSSIma(find(lSSIma(:,:,k))+read*fullphase*(k-1)));
end

muscleNSIma = zeros(read,fullphase,repeats,7);
muscleSSIma = zeros(read,fullphase,repeats,7);
muscleNSIma(:,:,:,1) = mgNSIma;
muscleSSIma(:,:,:,1) = mgSSIma;
muscleNSIma(:,:,:,2) = lgNSIma;
muscleSSIma(:,:,:,2) = lgSSIma;
muscleNSIma(:,:,:,3) = gNSIma;
muscleSSIma(:,:,:,3) = gSSIma;
muscleNSIma(:,:,:,4) = sNSIma;
muscleSSIma(:,:,:,4) = sSSIma;
muscleNSIma(:,:,:,5) = pNSIma;
muscleSSIma(:,:,:,5) = pSSIma;
muscleNSIma(:,:,:,6) = tNSIma;
muscleSSIma(:,:,:,6) = tSSIma;
muscleNSIma(:,:,:,7) = lNSIma;
muscleSSIma(:,:,:,7) = lSSIma;
% clear gNSIma gSSIma sNSIma sSSIma pNSIma pSSIma tNSIma tSSIma
Perf = zeros(repeats,1);
Perf_Img_tc = zeros(repeats,1);

for k = 2:repeats-1
    Perf(k,1) = -lambda/T*log(((SS(k,1)-(NS(k,1)+NS(k+1,1))/2)/(SS(k,1)+(NS(k,1)+NS(k+1,1))/2))*(1-exp(T/T1t))+1);
    Perf(k,2) = -lambda/T*log(((SS(k,2)-(NS(k,2)+NS(k+1,2))/2)/(SS(k,2)+(NS(k,2)+NS(k+1,2))/2))*(1-exp(T/T1t))+1);
    Perf(k,3) = -lambda/T*log(((SS(k,3)-(NS(k,3)+NS(k+1,3))/2)/(SS(k,3)+(NS(k,3)+NS(k+1,3))/2))*(1-exp(T/T1t))+1);
    Perf(k,4) = -lambda/T*log(((SS(k,4)-(NS(k,4)+NS(k+1,4))/2)/(SS(k,4)+(NS(k,4)+NS(k+1,4))/2))*(1-exp(T/T1t))+1);
    Perf(k,5) = -lambda/T*log(((SS(k,5)-(NS(k,5)+NS(k+1,5))/2)/(SS(k,5)+(NS(k,5)+NS(k+1,5))/2))*(1-exp(T/T1t))+1);
    Perf(k,6) = -lambda/T*log(((SS(k,6)-(NS(k,6)+NS(k+1,6))/2)/(SS(k,6)+(NS(k,6)+NS(k+1,6))/2))*(1-exp(T/T1t))+1);
    Perf(k,7) = -lambda/T*log(((SS(k,7)-(NS(k,7)+NS(k+1,7))/2)/(SS(k,7)+(NS(k,7)+NS(k+1,7))/2))*(1-exp(T/T1t))+1);
end

Perf2 = Perf;
% Perf(31:54,:) = NaN; % set exercise portion equal to NaN;

for k = 1:repeats
    for m = 1:7
        if Perf(k,m)>400
            Perf(k,m) = 0;
%             Perf(k,m)=(Perf(k-1,m)+Perf(k+1,m))/2;
        end
        if Perf(k,m)<-10
            Perf(k,m) = 0;
%             Perf(k,m)=(Perf(k-1,m)+Perf(k+1,m))/2;
        end
    end
end
%%
for k = 1:7
    Perf_smooth1(:,k) = smoothdata(Perf(:,k));
    Perf_smooth2(:,k) = smoothdata(Perf(:,k),'sgolay');
end

figure, plot(4:4:repeats*4,Perf)
figure, plot(4:4:repeats*4,Perf_smooth1)
figure, plot(4:4:repeats*4,Perf_smooth2)

vars.ASL.mask.mgMask = mgMask;
vars.ASL.mask.lgMask = lgMask;
vars.ASL.mask.gMask = gMask;
vars.ASL.mask.sMask = sMask;
vars.ASL.mask.pMask = pMask;
vars.ASL.mask.tMask = tMask;
vars.ASL.mask.lMask = lMask;

vars.ASL.Perf.AverageTimecourse = Perf;
vars.ASL.Perf.AverageTimecourse_smooth1 = Perf_smooth1;
vars.ASL.Perf.AverageTimecourse_smooth2 = Perf_smooth2;
%%
Perf_Img = zeros(size(ImaNS));
for iread = 1:read
    for iphase = 1:fullphase
        for irep = 2:repeats-1
            if sumM_threshold(iread,iphase) == 1;
                Perf_Img(iread,iphase,irep) = -lambda/T*log(((ImaSS(iread,iphase,irep)-(ImaNS(iread,iphase,irep)+ImaNS(iread,iphase,irep+1))/2)/(ImaSS(iread,iphase,irep)+(ImaNS(iread,iphase,irep)+ImaNS(iread,iphase,irep+1))/2))*(1-exp(T/T1t))+1);
                if Perf_Img(iread,iphase,irep)<0
                    Perf_Img(iread,iphase,irep) = NaN;
                elseif Perf_Img(iread,iphase,irep)>250
                    Perf_Img(iread,iphase,irep) = NaN;
                end
            end
        end
    end
end
for k = 1:repeats-1
    Perf_Img_temp = Perf_Img(:,:,k);
    temp2 = Perf_Img_temp(find(Perf_Img_temp));
    size_temp(k) = length(temp2)-sum(isnan(temp2));

    Perf_Img_tc(k,1) = mean(Perf_Img_temp(find(mgMask)),"omitnan");
    Perf_Img_tc(k,2) = mean(Perf_Img_temp(find(lgMask)),"omitnan");
    Perf_Img_tc(k,3) = mean(Perf_Img_temp(find(gMask)),"omitnan");
    Perf_Img_tc(k,4) = mean(Perf_Img_temp(find(sMask)),"omitnan");
    Perf_Img_tc(k,5) = mean(Perf_Img_temp(find(pMask)),"omitnan");
    Perf_Img_tc(k,6) = mean(Perf_Img_temp(find(tMask)),"omitnan");
    Perf_Img_tc(k,7) = mean(Perf_Img_temp(find(lMask)),"omitnan");


    Perf_Img_tc_med(k,1) = median(Perf_Img_temp(find(mgMask)),"omitnan");
    Perf_Img_tc_med(k,2) = median(Perf_Img_temp(find(lgMask)),"omitnan");
    Perf_Img_tc_med(k,3) = median(Perf_Img_temp(find(gMask)),"omitnan");
    Perf_Img_tc_med(k,4) = median(Perf_Img_temp(find(sMask)),"omitnan");
    Perf_Img_tc_med(k,5) = median(Perf_Img_temp(find(pMask)),"omitnan");
    Perf_Img_tc_med(k,6) = median(Perf_Img_temp(find(tMask)),"omitnan");
    Perf_Img_tc_med(k,7) = median(Perf_Img_temp(find(lMask)),"omitnan");
end
figure, plot(Perf_Img_tc), legend('medial gastroc','lateral gastroc','mean gastroc','soleus','peroneus','ta','leg')
figure, plot(Perf_Img_tc_med), legend('medial gastroc','lateral gastroc','mean gastroc','soleus','peroneus','ta','leg')
%%
figure, for k = 1:repeats-1, subplot(1,2,1), imagesc([ImaNS(:,:,k), ImaSS(:,:,k)]),axis image 
    subplot(1,2,2), imagesc(Perf_Img(:,:,k)), axis image, caxis([0 150]),colorbar, title(k), hold on, contour(gMask,'k'), contour(sMask,'k'), contour(pMask,'k'),contour(tMask,'k'), hold off, pause(0.1),drawnow, end
figure, plot(4:4:size(Perf,1)*4,Perf), legend('medial gastroc','lateral gastroc','mean gastroc','soleus','peroneus','ta','leg'), title('Perfusion Timecourse - ROI average')
figure, plot(4:4:size(Perf,1)*4,Perf_Img_tc), legend('medial gastroc','lateral gastroc','mean gastroc','soleus','peroneus','ta','leg'), title('Perfusion Timecourse - from perfusion image')
figure, plot(4:4:size(Perf,1)*4,Perf_Img_tc-Perf), legend('medial gastroc','lateral gastroc','mean gastroc','soleus','peroneus','ta','leg'), title('Perfusion Timecourse - from perfusion image')
% figure, plot(PerfImg_gAvg)

vars.ASL.Perf.Img = Perf_Img;

vars.ASL.Perf.AverageTimecourse_IMG_mean = Perf_Img_tc;
vars.ASL.Perf.AverageTimecourse_IMG_med = Perf_Img_tc_med;