function vPIVOT__Sup_PC_reorganization
% Reorganize Superior GRE data from both reference scans and keyhold
% acquisitions
%%
global vars

RefPre_raw = vars.PC.raw.Ref1;
RefPost_raw = vars.PC.raw.Ref2;
KHdata_raw = vars.PC.raw.KH;

channels = size(RefPre_raw,3);

RefPre = zeros(96,96,2,3,channels,2); % dimensions = read,phase,echoes,repeats,channels,VENC+/-
RefPost = zeros(96,96,2,5,channels,2);
KHdata = zeros(96,24,2,channels,vars.reps*2,2); % dimensions = read,phase,echoes,channels,repeats,VENC+/-

counter = 1;
for iRep = 1:3
    for iPE = 1:96
        for iEcho = 1:2
            for iVENC = 1:2
                RefPre(:,iPE,iEcho,iRep,:,iVENC) = RefPre_raw(:,counter,:);    
                RefPost(:,iPE,iEcho,iRep,:,iVENC) = RefPost_raw(:,counter,:); 
                counter = counter+1;
            end
        end
    end
end

 
counter = 1;
for iPE = 1:24
    for iEcho = 1:2
        for iVENC = 1:2
            KHdata(:,iPE,iEcho,:,:,iVENC) = KHdata_raw(:,counter,:,:);
            counter = counter+1;
        end
    end
end
%%
KHdata = permute(KHdata,[1 2 3 5 4 6]); % reorder as read, phase, echo, repeat, channel, VENC+/-
KHdata = KHdata(:,:,:,2:2:end,:,:); % discard data collected after NS inversion (odd acquisitions)

Data_wPreRef = zeros(96,96,2,vars.reps,channels,2);
Data_wPostRef = zeros(96,96,2,vars.reps,channels,2);
for k = 1:vars.reps
    Data_wPreRef(:,:,:,k,:,:) = mean(RefPre(:,:,:,2:3,:,:),4);
    Data_wPostRef(:,:,:,k,:,:) = mean(RefPost(:,:,:,2:3,:,:),4);
end
Data_wPreRef(:,37:60,:,:,:,:) = KHdata;
Data_wPostRef(:,37:60,:,:,:,:) = KHdata;

vars.PC.ksp.KHdata = KHdata;
vars.PC.ksp.PreRef = RefPre;
vars.PC.ksp.PostRef = RefPost;
vars.PC.ksp.Data_wPreRef = Data_wPreRef;
vars.PC.ksp.Data_wPostRef = Data_wPostRef;
