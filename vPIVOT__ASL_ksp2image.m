function vPIVOT__ASL_ksp2image
%%
global vars

read = 80;
phase = 50;
fullphase = 80;
repeats = vars.reps;
channels = size(vars.ASL.raw.EPI,3);

NSraw = vars.ASL.raw.EPI(:,:,:,1:2:end);
SSraw = vars.ASL.raw.EPI(:,:,:,2:2:end);

NSraw = permute(NSraw,[1 2 4 3]);
SSraw = permute(SSraw,[1 2 4 3]);

NSeven = zeros(size(NSraw));
SSeven = zeros(size(SSraw));
NSodd = zeros(size(NSraw));
SSodd = zeros(size(SSraw));

NSeven(:,2:2:phase,:,:) = NSraw([1 read:-1:2],2:2:phase,:,:);
SSeven(:,2:2:phase,:,:) = SSraw([1 read:-1:2],2:2:phase,:,:);
NSodd(:,1:2:phase-1,:,:) = NSraw(:,1:2:phase-1,:,:);
SSodd(:,1:2:phase-1,:,:) = SSraw(:,1:2:phase-1,:,:);   

KspPartNS = NSeven+NSodd;
KspPartSS = SSeven+SSodd;

clear NSraw NSeven NSodd SSraw SSeven SSodd

%Fill in k-space with complex conjugate
KspConjNS = conj(KspPartNS);
KspConjSS = conj(KspPartSS);
KspNS = zeros(read,fullphase,repeats,channels);
KspSS = zeros(read,fullphase,repeats,channels);

KspNS(:,1:fullphase-phase,:,:) = -KspConjNS(read:-1:1,phase:-1:2*phase-fullphase+1,:,:);
KspSS(:,1:fullphase-phase,:,:) = -KspConjSS(read:-1:1,phase:-1:2*phase-fullphase+1,:,:);

clear KspConjNS KspConjSS

KspNS(:,fullphase-phase+1:fullphase,:,:) = KspPartNS;
KspSS(:,fullphase-phase+1:fullphase,:,:) = KspPartSS;

clear KspPartNS KspPartSS

ImaNSComplex = fftshift(ifft2(ifftshift(KspNS)));
% write_mda('ImaNSComplex.mda',ImaNSComplex);
ImaNS = squeeze(sqrt(sum(abs(ImaNSComplex).^2,4)));
% write_mda('ImaNS.mda',ImaNS);
clear ImaNSComplex 

ImaSSComplex = fftshift(ifft2(ifftshift(KspSS)));
% write_mda('ImaSSComplex.mda',ImaSSComplex);
ImaSS = squeeze(sqrt(sum(abs(ImaSSComplex).^2,4)));
% write_mda('ImaSS.mda',ImaSS);
clear ImaSSComplex 

%


figure(1), subplot(2,2,1), imagesc(ImaNS(:,:,10)), axis image, colormap gray, title('NS before correction')
subplot(2,2,2), imagesc(ImaSS(:,:,10)), axis image, colormap gray, title('SS before correction')
% 

refraw = vars.ASL.raw.refEPI; % expectation is dimensions is read,phase,repeats,channels
refraw = permute(refraw,[1 2 4 3]);

% Perform N/2 ghost correction
refeven = squeeze(refraw([1 read:-1:2],2:2:phase,1,:));
refodd = squeeze(refraw(:,1:2:phase-1,1,:));
Imarefeven = fftshift(ifft(ifftshift(refeven)));
Imarefodd = fftshift(ifft(ifftshift(refodd)));
EPIphase = squeeze(angle(sum(Imarefodd.*conj(Imarefeven),3)));
for k = 1:80
    for j = 1:25
        if EPIphase(k,j) < -2
            EPIphase(k,j) = EPIphase(k,j)+2*pi;
        end
    end
end

% clear refeven refodd Imarefeven Imarefodd refraw
figure, imagesc(EPIphase), impixelinfo
%
[bw,xref,yref]= roipoly;
phasecorrection = zeros(3,phase/2);
for k = 1:phase/2
    phasecorrection(:,k) = polyfit(ceil(min(yref)):floor(max(yref)),EPIphase(ceil(min(yref)):floor(max(yref)),k)',2);
end
phasecorrectionavg = mean(phasecorrection,2);
clear phasecorrection
figure, hold on
plot(ceil(min(yref)):floor(max(yref)), EPIphase(ceil(min(yref)):floor(max(yref)),:)')
plot(ceil(min(yref)):floor(max(yref)), (ceil(min(yref)):floor(max(yref))).^2*phasecorrectionavg(1)+(ceil(min(yref)):floor(max(yref)))*phasecorrectionavg(2)+phasecorrectionavg(3),'.')
correctphase = zeros(read,fullphase,repeats,channels);

for k = 1:read
     correctphase(k,1:2:fullphase-1,:,:) = k^2*phasecorrectionavg(1)+k*phasecorrectionavg(2)+phasecorrectionavg(3);
end
clear phasecorrectionavg

NSodd_correctphase = zeros(read,fullphase,repeats,channels);
NSeven_correctphase = zeros(read,fullphase,repeats,channels);
NSodd_correctphase(:,1:2:end,:,:) = KspNS(:,1:2:end,:,:);
NSeven_correctphase(:,2:2:end,:,:) = KspNS(:,2:2:end,:,:);

proj_NSodd = fftshift(ifft(ifftshift(NSodd_correctphase)));
clear NSodd_correctphase
NSodd_correct = exp(-1i*correctphase).*proj_NSodd;
clear proj_NSodd
Ksp_NSodd_correct = ifftshift(fft(fftshift(NSodd_correct)));
clear NSodd_correct
KspNS = Ksp_NSodd_correct+NSeven_correctphase;
clear NSeven_correctphase Ksp_NSodd_correct

ImaNScomplex = zeros(size(KspNS));
for k = 1:repeats
    ImaNScomplex(:,:,k,:) = fftshift(ifft2(ifftshift(KspNS(:,:,k,:))));
    ImaNS(:,:,k) = squeeze(sqrt(sum(abs(ImaNScomplex(:,:,k,:)).^2,4)));
end

% write_mda('ImaNS_NoMotionCorrection.mda',ImaNS);
% write_mda('ImaNScomplex.mda',ImaNScomplex);
clear ImaNScomplex KspNS

SSodd_correctphase = zeros(read,fullphase,repeats,channels);
SSeven_correctphase = zeros(read,fullphase,repeats,channels);
SSodd_correctphase(:,1:2:end,:,:) = KspSS(:,1:2:end,:,:);
SSeven_correctphase(:,2:2:end,:,:) = KspSS(:,2:2:end,:,:);

proj_SSodd = fftshift(ifft(ifftshift(SSodd_correctphase)));
clear SSodd_correctphase
SSodd_correct = exp(-1i*correctphase).*proj_SSodd;
clear proj_SSodd
Ksp_SSodd_correct = ifftshift(fft(fftshift(SSodd_correct)));
clear SSodd_correct
KspSS = Ksp_SSodd_correct+SSeven_correctphase;
clear SSeven_correctphase Ksp_SSodd_correct

ImaSScomplex = zeros(size(KspSS));
for k = 1:repeats
    ImaSScomplex(:,:,k,:) = fftshift(ifft2(ifftshift(KspSS(:,:,k,:))));
    ImaSS(:,:,k) = squeeze(sqrt(sum(abs(ImaSScomplex(:,:,k,:)).^2,4)));
end

% write_mda('ImaSS_NoMotionCorrection.mda',ImaSS);
% write_mda('ImaSScomplex.mda',ImaSScomplex);
clear ImaSScomplex KspSS
clear xref yref correctphase bw EPIphase


figure(1), subplot(2,2,3), imagesc(ImaNS(:,:,10)), axis image, colormap gray, title('NS after correction')
subplot(2,2,4), imagesc(ImaSS(:,:,10)), axis image, colormap gray, title('SS after correction')

vars.ASL.img.ImaNS = ImaNS;
vars.ASL.img.ImaSS = ImaSS;


%%
% UNCOMMENT TO PERFORM MOTION CORRECTION

% javaaddpath '/Users/erinenglund/Documents/MATLAB/ij.jar';
% javaaddpath '/Users/erinenglund/Documents/MATLAB/mij.jar';
% MIJ.start('/Users/erinenglund/Documents/ImageJ');
% 
% %%% -Check data for motion artifacts and Do Registeration if required.
% clear UnRegCombData
% %     UnRegCombData = cat(3, Img, Img20_Off);
% %      UnRegCombData = cat(3, Img1, Img2, Img20_Off, ImgFreqWassr, ImgB1);
% 
% ImaStack = zeros(read,fullphase,repeats*2);
% ImaStack(:,:,1:2:repeats*2) = ImaNS;
% ImaStack(:,:,2:2:repeats*2) = ImaSS;
% 
% UnRegCombData = ImaStack;
% for n = 1:size(UnRegCombData,3)
%     UnRegCombData1(:,:,n) = double(UnRegCombData(:,:,n));%.*(double(mask==2));
% end
% MIJ.createImage(UnRegCombData1);  % after this use 'stackreg' plugin
% 
% uiwait(msgbox('Use-''ImageJ->Plugins->stackreg->Rigid Body'' to perform Registeration'));
% 
% prompt={'Registeration status-> 0 (No) and 1 (Yes)'}; defans={'1'}; fields = {'num'}; options.Resize='on';
% info = inputdlg(prompt, 'Registeration status?', 1, defans, options);
% if ~isempty(info)              %see if user hit cancel
%     info = cell2struct(info,fields);   RegStatus = str2num(info.num);   %convert string to number
% else  msgbox('Cancelled??', 'Hey!')
% end
% if(RegStatus==0), msgbox('No Registeration done!!', 'Hey!'); end
% 
% if RegStatus>0, RegCombData = MIJ.getImage('Import from Matlab');
% else    RegCombData = UnRegCombData;
% end
% 
% MIJ.exit;
% 
% clear Img Img20_Off
% %     Img = RegCombData(:,:,1:nuPoint);
% %     Img20_Off = RegCombData(:,:,nuPoint+1:(nuPoint+length(filename20)));
% 
% ImgStackCorr = RegCombData;
% %     Img20_Off = RegCombData(:,:,nuPoint+1:nuPoint+2);
% %     ImgFreqWassr = RegCombData(:,:,nuPoint+1:nuPoint+22); 
% %     ImgB1 = RegCombData(:,:,nuPoint+23:nuPoint+24);  
% % end
% % 
% % ImaNS_Uncorr = ImaNS;
% % ImaSS_Uncorr = ImaSS;
% 
% ImaNS_MC = ImgStackCorr(:,:,1:2:repeats*2);
% ImaSS_MC = ImgStackCorr(:,:,2:2:repeats*2);
% 
% % write_mda('ImaSS.mda',ImaSS);
% % write_mda('ImaNS.mda',ImaNS);
% % 
% % figure
% % for k = 1:10:repeats
% %     subplot(1,2,1), imagesc(ImaNS(:,:,k)-ImaNS_Uncorr(:,:,k)), axis image, title('NS Corrected-Uncorrected')
% %     subplot(1,2,2), imagesc(ImaSS(:,:,k)-ImaSS_Uncorr(:,:,k)), axis image, title('SS Corrected-Uncorrected')
% %     drawnow, pause(0.1)
% % end
% 
% clear ImaNS_Uncorr ImaSS_Uncorr ImaStack ImgStackCorr RegCombData RegStatus UnRegCombData UnRegCombData1 defans fields info n options prompt

