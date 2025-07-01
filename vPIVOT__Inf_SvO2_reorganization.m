function vPIVOT__Inf_SvO2_reorganization
% Reorganize Inferior GRE data from both reference scans and keyhold
% acquisitions
%%
global vars

RefPre_raw = vars.InfGRE.raw.Ref1;
RefPost_raw = vars.InfGRE.raw.Ref2;
KHdata_raw = vars.InfGRE.raw.KH;

channels = size(RefPre_raw,3);
RefPre = zeros(96,96,5,3,channels);
RefPost = RefPre;
KHdata = zeros(96,24,5,channels,vars.reps*2);

counter = 1;
for iRep = 1:3
    for iPE = 1:96
        for iEcho = 1:5
            RefPre(:,iPE,iEcho,iRep,:) = RefPre_raw(:,counter,:);
            RefPost(:,iPE,iEcho,iRep,:) = RefPost_raw(:,counter,:);
            counter = counter+1;
        end
    end
end

counter = 1;
for iPE = 1:24
    for iEcho = 1:5
        KHdata(:,iPE,iEcho,:,:) = KHdata_raw(:,counter,:,:);
        counter = counter+1;
    end
end
KHdata = permute(KHdata,[1 2 3 5 4]); % reorder as read, phase, echo, repeat, channel
KHdata = KHdata(:,:,:,2:2:end,:); % discard data collected after NS inversion (odd acquisitions)

Data_wPreRef = zeros(96,96,5,vars.reps,channels);
Data_wPostRef = zeros(96,96,5,vars.reps,channels);
for k = 1:vars.reps
    Data_wPreRef(:,:,:,k,:) = mean(RefPre(:,:,:,2:3,:),4);
    Data_wPostRef(:,:,:,k,:) = mean(RefPost(:,:,:,2:3,:),4);
end
%%
Data_wPreRef(:,38:60,:,:,:) = KHdata(:,2:end,:,:,:);
Data_wPostRef(:,38:60,:,:,:) = KHdata(:,2:end,:,:,:);
%%
vars.InfGRE.ksp.KHdata = KHdata;
vars.InfGRE.ksp.PreRef = RefPre;
vars.InfGRE.ksp.PostRef = RefPost;
vars.InfGRE.ksp.Data_wPreRef = Data_wPreRef;
vars.InfGRE.ksp.Data_wPostRef = Data_wPostRef;
