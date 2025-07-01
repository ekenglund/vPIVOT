function vPIVOT__Read_vPIVOT_datXA30A(datfile,reps)

% University of Pennsylvania
% Written by Erin Englund
% March  2017 - Modified by Ana E. Rodríguez-Soto to read VE11C
% August 2017 - Modified by Ana E. Rodríguez-Soto to read OxBOLD  data (VE11C)
% January 2018 - Modified by Ana E. Rodríguez-Soto to read OxBOLD GRE and PC data
% March 2021 - Modified by Erin Englund to read vPIVOT data
% November 2022 - Modified by Erin Englund to read vPIVOT data XA30A

% Open *dat file and skip the header
global vars;
reps = reps*2;


%%
% Open *dat file and skip the header
fid = fopen(datfile,'r');
fread(fid,1,'uint32');
nMeas=fread(fid,1,'uint32');
measLen=0;
for iMeas=1:nMeas-1;
    fseek(fid,16,0);
    tmp=fread(fid,1,'uint64');
    measLen=measLen+ceil(tmp/512)*512;
    fseek(fid,152-24,0);
end
fseek(fid,2*4+152*64+126*4+measLen,-1);
headerSize = fread(fid,1,'uint32');
fread(fid,headerSize-4,'uint8');
fread(fid,48,'uint8');

% Find number of readouts and channels used.
num_readouts = fread(fid,1,'uint16');
num_channels = fread(fid,1,'uint16');
fseek(fid, -48-4, 0);

X = sprintf('Reading:   %s      \nProgress...',datfile);
disp(X);
tic;
percentFinished=0;
%
nviews_SUP_Ref1 = 96*2*2*3;  % = # PE lines (96) x VENCs x echoes x repeats
nviews_INF_Ref1 = 96*5*3;    % = # PE lines (96) x echoes x repeats
nviews_M0 = 50;              % = # PE lines for EPI M0 acquisition
% Looped Reps number of times
nviews_HbO2_KH = 24*5;        % = # PE lines (24) x echoes
nviews_EPI = 50;             % = # PE lines (50)
nviews_KHdummy = 2*2;        % = # 1 dummy for each VENC x echoes
nviews_PC_KH = 24*2*2;       % = # PE lines (24) x VENCs x echoes
% then only once
nviews_EPI_Ref = 50*2;       % = # PE lines (50) x NS and SS
nviews_SUP_Ref2 = 96*2*2*3;  % = # PE lines (96) x VENCs x echoes x repeats
nviews_INF_Ref2 = 96*5*3;    % = # PE lines (96) x echoes x repeats


SUP_Ref1 = zeros(96,nviews_SUP_Ref1,num_channels);
INF_Ref1 = zeros(96,nviews_INF_Ref1,num_channels);
M0 = zeros(80,nviews_M0,num_channels);

HbO2_KH = zeros(96,nviews_HbO2_KH,num_channels,reps);
EPI = zeros(80,nviews_EPI,num_channels,reps);
PC_KH_dummy = zeros(96,nviews_KHdummy,num_channels,reps);
PC_KH = zeros(96,nviews_PC_KH,num_channels,reps);

EPI_Ref = zeros(80,nviews_EPI_Ref,num_channels);
SUP_Ref2 = zeros(96,nviews_SUP_Ref2,num_channels);
INF_Ref2 = zeros(96,nviews_INF_Ref2,num_channels);

%% Read data.
% data is stored as an array of:
% (size of kx) x (#of total lines acquired, including all images and repeats) x (#of channels)
iPE = 0;
while 1
    iPE = iPE+1;
    
    fread(fid,48,'uint8');
    num_readouts=fread(fid,1,'uint16');
    
    if num_readouts == 0
        fseek(fid,-50,'cof');
        
        ulDMALength = 184;
        mdhStart = 1-ulDMALength;
        data_u8 = fread( fid, ulDMALength, 'uint8=>uint8' );
        data_u8 = data_u8( mdhStart+end :  end );
        
        data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
        ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );
        fread( fid, ulDMALength-184, 'uint8=>uint8' );
        iPE = iPE-1;
        continue;
    end
    
    fread(fid,192-50,'uint8');
    
    for iCh = 1:num_channels
        fread(fid,32,'uint8');
        tmp = fread(fid,num_readouts*2,'float32');
        tmp = tmp(1:2:end)+1i*tmp(2:2:end);
        SUP_Ref1(:,iPE,iCh)=tmp;
    end
    
    if (iPE == nviews_SUP_Ref1)
        disp('Superior Ref1 Done');
        break;
    end
end

iPE = 0;
while 1
    iPE = iPE+1;
    
    fread(fid,48,'uint8');
    num_readouts=fread(fid,1,'uint16');
    
    if num_readouts == 0
        fseek(fid,-50,'cof');
        
        ulDMALength = 184;
        mdhStart = 1-ulDMALength;
        data_u8 = fread( fid, ulDMALength, 'uint8=>uint8' );
        data_u8 = data_u8( mdhStart+end :  end );
        
        data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
        ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );
        fread( fid, ulDMALength-184, 'uint8=>uint8' );
        iPE = iPE-1;
        continue;
    end
    
    fread(fid,192-50,'uint8');
    
    for iCh = 1:num_channels
        fread(fid,32,'uint8');
        tmp = fread(fid,num_readouts*2,'float32');
        tmp = tmp(1:2:end)+1i*tmp(2:2:end);
        INF_Ref1(:,iPE,iCh)=tmp;
    end
    
    if (iPE == nviews_INF_Ref1)
        disp('Inferior Ref1 Done')
        break;
    end
end

iPE = 0;
while 1
    iPE = iPE+1;
    
    fread(fid,48,'uint8');
    num_readouts=fread(fid,1,'uint16');
    
    if num_readouts == 0
        fseek(fid,-50,'cof');
        
        ulDMALength = 184;
        mdhStart = 1-ulDMALength;
        data_u8 = fread( fid, ulDMALength, 'uint8=>uint8' );
        data_u8 = data_u8( mdhStart+end :  end );
        
        data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
        ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );
        fread( fid, ulDMALength-184, 'uint8=>uint8' );
        iPE = iPE-1;
        continue;
    end
    
    fread(fid,192-50,'uint8');
    
    for iCh = 1:num_channels
        fread(fid,32,'uint8');
        tmp = fread(fid,num_readouts*2,'float32');
        tmp = tmp(1:2:end)+1i*tmp(2:2:end);
        M0(:,iPE,iCh)=tmp;
    end
    
    if (iPE == nviews_M0)
        disp('M0 Done')
        break;
    end
end

disp('Dynamic vPIVOT...')

for iRep = 1:reps
    iPE = 0;
    while 1
        iPE = iPE+1;
        
        fread(fid,48,'uint8');
        num_readouts=fread(fid,1,'uint16');
        
        if num_readouts == 0
            fseek(fid,-50,'cof');
            
            ulDMALength = 184;
            mdhStart = 1-ulDMALength;
            data_u8 = fread( fid, ulDMALength, 'uint8=>uint8' );
            data_u8 = data_u8( mdhStart+end :  end );
            
            data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
            ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );
            fread( fid, ulDMALength-184, 'uint8=>uint8' );
            iPE = iPE-1;
            continue;
        end
        
        fread(fid,192-50,'uint8');
        
        for iCh = 1:num_channels
            fread(fid,32,'uint8');
            tmp = fread(fid,num_readouts*2,'float32');
            tmp = tmp(1:2:end)+1i*tmp(2:2:end);
            HbO2_KH(:,iPE,iCh,iRep)=tmp;
        end
        
        if (iPE == nviews_HbO2_KH)
            break;
        end
    end

    iPE = 0;
    while 1
        iPE = iPE+1;
        
        fread(fid,48,'uint8');
        num_readouts=fread(fid,1,'uint16');
        
        if num_readouts == 0
            fseek(fid,-50,'cof');
            
            ulDMALength = 184;
            mdhStart = 1-ulDMALength;
            data_u8 = fread( fid, ulDMALength, 'uint8=>uint8' );
            data_u8 = data_u8( mdhStart+end :  end );
            
            data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
            ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );
            fread( fid, ulDMALength-184, 'uint8=>uint8' );
            iPE = iPE-1;
            continue;
        end
        
        fread(fid,192-50,'uint8');
        
        for iCh = 1:num_channels
            fread(fid,32,'uint8');
            tmp = fread(fid,num_readouts*2,'float32');
            tmp = tmp(1:2:end)+1i*tmp(2:2:end);
            EPI(:,iPE,iCh,iRep)=tmp;
        end
        
        if (iPE == nviews_EPI)
            break;
        end
    end
    
    iPE = 0;
    while 1
        iPE = iPE+1;
        
        fread(fid,48,'uint8');
        num_readouts=fread(fid,1,'uint16');
        
        if num_readouts == 0
            fseek(fid,-50,'cof');
            
            ulDMALength = 184;
            mdhStart = 1-ulDMALength;
            data_u8 = fread( fid, ulDMALength, 'uint8=>uint8' );
            data_u8 = data_u8( mdhStart+end :  end );
            
            data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
            ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );
            fread( fid, ulDMALength-184, 'uint8=>uint8' );
            iPE = iPE-1;
            continue;
        end
        
        fread(fid,192-50,'uint8');
        
        for iCh = 1:num_channels
            fread(fid,32,'uint8');
            tmp = fread(fid,num_readouts*2,'float32');
            tmp = tmp(1:2:end)+1i*tmp(2:2:end);
            PC_KH_dummy(:,iPE,iCh,iRep)=tmp;
        end
        
        if (iPE == nviews_KHdummy)
            break;
        end
    end
    
    iPE = 0;
    while 1
        iPE = iPE+1;
        
        fread(fid,48,'uint8');
        num_readouts=fread(fid,1,'uint16');
        
        if num_readouts == 0
            fseek(fid,-50,'cof');
            
            ulDMALength = 184;
            mdhStart = 1-ulDMALength;
            data_u8 = fread( fid, ulDMALength, 'uint8=>uint8' );
            data_u8 = data_u8( mdhStart+end :  end );
            
            data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
            ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );
            fread( fid, ulDMALength-184, 'uint8=>uint8' );
            iPE = iPE-1;
            continue;
        end
        
        fread(fid,192-50,'uint8');
        
        for iCh = 1:num_channels
            fread(fid,32,'uint8');
            tmp = fread(fid,num_readouts*2,'float32');
            tmp = tmp(1:2:end)+1i*tmp(2:2:end);
            PC_KH(:,iPE,iCh,iRep)=tmp;
        end
        
        if (iPE == nviews_PC_KH)
            break;
        end
    end
    if floor((100*iRep)/reps) > percentFinished+9
        percentFinished = floor((100*iRep)/reps);
        progress_str    = sprintf('  %3.0f %%  ', percentFinished);
        fprintf('%s \n',progress_str);
    end
end
disp('Dynamic vPIVOT Done')
    
iPE = 0;
while 1
    iPE = iPE+1;
    
    fread(fid,48,'uint8');
    num_readouts=fread(fid,1,'uint16');
    
    if num_readouts == 0
        fseek(fid,-50,'cof');
        
        ulDMALength = 184;
        mdhStart = 1-ulDMALength;
        data_u8 = fread( fid, ulDMALength, 'uint8=>uint8' );
        data_u8 = data_u8( mdhStart+end :  end );
        
        data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
        ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );
        fread( fid, ulDMALength-184, 'uint8=>uint8' );
        iPE = iPE-1;
        continue;
    end
    
    fread(fid,192-50,'uint8');
    
    for iCh = 1:num_channels
        fread(fid,32,'uint8');
        tmp = fread(fid,num_readouts*2,'float32');
        tmp = tmp(1:2:end)+1i*tmp(2:2:end);
        EPI_Ref(:,iPE,iCh)=tmp;
    end
    
    if (iPE == nviews_EPI_Ref)
        disp('EPI Ref Done');
        break;
    end
end
    
iPE = 0;
while 1
    iPE = iPE+1;
    
    fread(fid,48,'uint8');
    num_readouts=fread(fid,1,'uint16');
    
    if num_readouts == 0
        fseek(fid,-50,'cof');
        
        ulDMALength = 184;
        mdhStart = 1-ulDMALength;
        data_u8 = fread( fid, ulDMALength, 'uint8=>uint8' );
        data_u8 = data_u8( mdhStart+end :  end );
        
        data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
        ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );
        fread( fid, ulDMALength-184, 'uint8=>uint8' );
        iPE = iPE-1;
        continue;
    end
    
    fread(fid,192-50,'uint8');
    
    for iCh = 1:num_channels
        fread(fid,32,'uint8');
        tmp = fread(fid,num_readouts*2,'float32');
        tmp = tmp(1:2:end)+1i*tmp(2:2:end);
        SUP_Ref2(:,iPE,iCh)=tmp;
    end
    
    if (iPE == nviews_SUP_Ref2)
        disp('Superior Ref2 Done');
        break;
    end
end
    
iPE = 0;
while 1
    iPE = iPE+1;
    
    fread(fid,48,'uint8');
    num_readouts=fread(fid,1,'uint16');
    
    if num_readouts == 0
        fseek(fid,-50,'cof');
        
        ulDMALength = 184;
        mdhStart = 1-ulDMALength;
        data_u8 = fread( fid, ulDMALength, 'uint8=>uint8' );
        data_u8 = data_u8( mdhStart+end :  end );
        
        data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
        ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );
        fread( fid, ulDMALength-184, 'uint8=>uint8' );
        iPE = iPE-1;
        continue;
    end
    
    fread(fid,192-50,'uint8');
    
    for iCh = 1:num_channels
        fread(fid,32,'uint8');
        tmp = fread(fid,num_readouts*2,'float32');
        tmp = tmp(1:2:end)+1i*tmp(2:2:end);
        INF_Ref2(:,iPE,iCh)=tmp;
    end
    
    if (iPE == nviews_INF_Ref2)
        disp('Inferior Ref2 Done');
        break;
    end
end




toc
fclose(fid);


    
    
    
    
    
%     
%     %%
%     SUP_Ref1 = zeros(num_readouts,96*2*2*3,num_channels);
%     for iPE = 1:96*2*2*3 % = # PE lines (96) x VENCs x echoes x repeats
%         fread(fid,48,'uint8');
%         num_readouts=fread(fid,1,'uint16');
%         
%         if num_readouts == 0
%             fseek(fid,-50,'cof');
%             
%             ulDMALength = 184;
%             mdhStart = 1-ulDMALength;
%             data_u8 = fread( fid, ulDMALength, 'uint8=>uint8' );
%             data_u8 = data_u8( mdhStart+end :  end );
%             
%             data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
%             ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );
%             fread( fid, ulDMALength-184, 'uint8=>uint8' );
%             iPE = iPE-1;
%             continue;
%         end
%         
%         fread(fid,192-50,'uint8');
%         
%         for iCh = 1:num_channels
%             fread(fid,32,'uint8');
%             tmp = fread(fid,num_readouts*2,'float32');
%             tmp = tmp(1:2:end)+1i*tmp(2:2:end);
%             SUP_Ref1(:,iPE,iCh)=tmp;
%         end
%     end
%     
%     
%     for iPE = 1:96*5*3 % = # PE lines (96) x echoes x repeats
%         fread(fid,48,'uint8');
%         num_readouts=fread(fid,1,'uint16');
%         if num_readouts == 0
%             fseek(fid,-50,'cof');
%             
%             ulDMALength = 184;
%             mdhStart = 1-ulDMALength;
%             data_u8 = fread( fid, ulDMALength, 'uint8=>uint8' );
%             data_u8 = data_u8( mdhStart+end :  end );
%             
%             data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
%             ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );
%             fread( fid, ulDMALength-184, 'uint8=>uint8' );
%             iPE = iPE-1;
%             continue;
%         end
%         
%         fread(fid,192-50,'uint8');
%         for iCh = 1:num_channels
%             fread(fid,32,'uint8');
%             tmp = fread(fid,num_readouts*2,'float32');
%             tmp = tmp(1:2:end)+1i*tmp(2:2:end);
%             INF_Ref1(:,iPE,iCh)=tmp;
%             clear tmp
%         end
%     end
%     %%
%     % ASL M0
%     for iPE = 1:50 % 50 = # PE lines for EPI M0 acquisition
%         fread(fid,48,'uint8');
%         num_readouts=fread(fid,1,'uint16');
%         
%         if num_readouts == 0
%             fseek(fid,-50,'cof');
%             
%             ulDMALength = 184;
%             mdhStart = 1-ulDMALength;
%             data_u8 = fread( fid, ulDMALength, 'uint8=>uint8' );
%             data_u8 = data_u8( mdhStart+end :  end );
%             
%             data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
%             ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );
%             fread( fid, ulDMALength-184, 'uint8=>uint8' );
%             iPE = iPE-1;
%             continue;
%         end
%         
%         fread(fid,192-50,'uint8');
%         for iCh = 1:num_channels
%             fread(fid,32,'uint8');
%             tmp = fread(fid,num_readouts*2,'float32');
%             tmp = tmp(1:2:end)+1i*tmp(2:2:end);
%             M0(:,iPE,iCh)=tmp;
%         end
%     end
    %%
    %%
%     
%     
%     data = zeros(num_readouts,nviews,num_channels);
%     %% Read data.
%     % data is stored as an array of:
%     % (size of kx) x (#of total lines acquired, including all images and repeats) x (#of channels)
%     iPE = 0;
%     while 1
%         iPE = iPE+1;
%         
%         
%         fread(fid,48,'uint8');
%         num_readouts=fread(fid,1,'uint16');
%         
%         if num_readouts == 0
%             fseek(fid,-50,'cof');
%             
%             ulDMALength = 184;
%             mdhStart = 1-ulDMALength;
%             data_u8 = fread( fid, ulDMALength, 'uint8=>uint8' );
%             data_u8 = data_u8( mdhStart+end :  end );
%             
%             data_u8(4)= bitget( data_u8(4),1);  % ubit24: keep only 1 bit from the 4th byte
%             ulDMALength = double( typecast( data_u8(1:4), 'uint32' ) );
%             fread( fid, ulDMALength-184, 'uint8=>uint8' );
%             iPE = iPE-1;
%             continue;
%         end
%         
%         
%         fread(fid,192-50,'uint8');
%         
%         for iCh = 1:num_channels
%             fread(fid,32,'uint8');
%             tmp = fread(fid,num_readouts*2,'float32');
%             tmp = tmp(1:2:end)+1i*tmp(2:2:end);
%             data(:,iPE,iCh)=tmp;
%         end
%         
%         if floor((100*iPE)/nviews) > percentFinished+9
%             percentFinished = floor((100*iPE)/nviews);
%             progress_str    = sprintf('  %3.0f %%  ', percentFinished);
%             fprintf('%s \n',progress_str);
%         end
%         
%         if (iPE == nviews)
%             break;
%         end
%     end
%     toc
%     fclose(fid);
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     % Read M0 data.
%     % data is stored as an array of:
%     % (size of kx) x (#of total lines acquired, including all images and repeats) x (#of channels)
%     % PIVOT acquisition structure:
%     % M0 (1 repeat, 80x50 matrix)
%     % HbO2 Base (5 repeats, 5 echoes, 96x96 matrix)
%     % HbO2 Ref1 (3 repeats, 5 echoes, 96x96 matrix)
%     % PIVOT (user set repeats, KH multi-echo GRE (5 echoes, 96x24 matrix)
%     %           --> EPI readout (80x50 matrix))
%     % Phase correction EPI readout (2 repeats, 80x50 matrix)
%     % HbO2 Ref2 (3 repeats, 5 echoes, 96x96 matrix)
%     
%     
%     % REF1 SUP PC dual echo GRE Baseline
%     for iPE = 1:96*2*2*3 % = # PE lines (96) x VENCs x echoes x repeats
%         fread(fid,48,'uint8');
%         num_readouts=fread(fid,1,'uint16');
%         fread(fid,192-50,'uint8');
%         for iCh = 1:num_channels
%             fread(fid,32,'uint8');
%             tmp = fread(fid,num_readouts*2,'float32');
%             tmp = tmp(1:2:end)+1i*tmp(2:2:end);
%             SUP_Ref1(:,iPE,iCh)=tmp;
%             clear tmp
%         end
%     end
%     %
%     % REF1 INF multi-echo GRE Baseline
%     for iPE = 1:96*5*3 % = # PE lines (96) x echoes x repeats
%         fread(fid,48,'uint8');
%         num_readouts=fread(fid,1,'uint16');
%         fread(fid,192-50,'uint8');
%         for iCh = 1:num_channels
%             fread(fid,32,'uint8');
%             tmp = fread(fid,num_readouts*2,'float32');
%             tmp = tmp(1:2:end)+1i*tmp(2:2:end);
%             INF_Ref1(:,iPE,iCh)=tmp;
%             clear tmp
%         end
%     end
%     %
%     % ASL M0
%     for iPE = 1:50 % 50 = # PE lines for EPI M0 acquisition
%         fread(fid,48,'uint8');
%         num_readouts=fread(fid,1,'uint16');
%         fread(fid,192-50,'uint8');
%         for iCh = 1:num_channels
%             fread(fid,32,'uint8');
%             tmp = fread(fid,num_readouts*2,'float32');
%             tmp = tmp(1:2:end)+1i*tmp(2:2:end);
%             M0(:,iPE,iCh)=tmp;
%         end
%     end
%     %%
%     
%     % PIVOT (KH Multi-echo GRE, then EPI, then dual-echo PC
%     % 2X number of repeats since 1 scan repeat comprises both NS and SS preps
%     
%     for iRep = 1:reps
%         for iPE = 1:24*5 % = # PE lines (24) x echoes
%             fread(fid,48,'uint8');
%             num_readouts=fread(fid,1,'uint16');
%             fread(fid,192-50,'uint8');
%             for iCh = 1:num_channels
%                 fread(fid,32,'uint8');
%                 tmp = fread(fid,num_readouts*2,'float32');
%                 tmp = tmp(1:2:end)+1i*tmp(2:2:end);
%                 HbO2_KH(:,iPE,iCh,iRep)=tmp;
%                 clear tmp
%             end
%         end
%         for iPE = 1:50 % = # PE lines (50)
%             fread(fid,48,'uint8');
%             num_readouts=fread(fid,1,'uint16');
%             fread(fid,192-50,'uint8');
%             for iCh = 1:num_channels
%                 fread(fid,32,'uint8');
%                 tmp = fread(fid,num_readouts*2,'float32');
%                 tmp = tmp(1:2:end)+1i*tmp(2:2:end);
%                 EPI(:,iPE,iCh,iRep)=tmp;
%                 clear tmp
%             end
%         end
%         
%         for iPE = 1:2*2 % = # 1 dummy for each VENC x echoes
%             fread(fid,48,'uint8');
%             num_readouts=fread(fid,1,'uint16');
%             fread(fid,192-50,'uint8');
%             for iCh = 1:num_channels
%                 fread(fid,32,'uint8');
%                 tmp = fread(fid,num_readouts*2,'float32');
%                 tmp = tmp(1:2:end)+1i*tmp(2:2:end);
%                 PC_KH_dummy(:,iPE,iCh,iRep)=tmp;
%                 clear tmp
%             end
%         end
%         
%         for iPE = 1:24*2*2 % = # PE lines (24) x VENCs x echoes
%             fread(fid,48,'uint8');
%             num_readouts=fread(fid,1,'uint16');
%             fread(fid,192-50,'uint8');
%             for iCh = 1:num_channels
%                 fread(fid,32,'uint8');
%                 tmp = fread(fid,num_readouts*2,'float32');
%                 tmp = tmp(1:2:end)+1i*tmp(2:2:end);
%                 PC_KH(:,iPE,iCh,iRep)=tmp;
%                 clear tmp
%             end
%         end
%     end
%     %%
%     
%     for iPE = 1:50*2 % = # PE lines (50) x NS and SS
%         fread(fid,48,'uint8');
%         num_readouts=fread(fid,1,'uint16');
%         fread(fid,192-50,'uint8');
%         for iCh = 1:num_channels
%             fread(fid,32,'uint8');
%             tmp = fread(fid,num_readouts*2,'float32');
%             tmp = tmp(1:2:end)+1i*tmp(2:2:end);
%             EPI_Ref(:,iPE,iCh)=tmp;
%             clear tmp
%         end
%     end
%     %%
%     % REF1 SUP PC dual echo GRE Baseline
%     for iPE = 1:96*2*2*3 % = # PE lines (96) x VENCs x echoes x repeats
%         fread(fid,48,'uint8');
%         num_readouts=fread(fid,1,'uint16');
%         fread(fid,192-50,'uint8');
%         for iCh = 1:num_channels
%             fread(fid,32,'uint8');
%             tmp = fread(fid,num_readouts*2,'float32');
%             tmp = tmp(1:2:end)+1i*tmp(2:2:end);
%             SUP_Ref2(:,iPE,iCh)=tmp;
%             clear tmp
%         end
%     end
%     %
%     % REF1 INF multi-echo GRE Baseline
%     for iPE = 1:96*5*3 % = # PE lines (96) x echoes x repeats
%         fread(fid,48,'uint8');
%         num_readouts=fread(fid,1,'uint16');
%         fread(fid,192-50,'uint8');
%         for iCh = 1:num_channels
%             fread(fid,32,'uint8');
%             tmp = fread(fid,num_readouts*2,'float32');
%             tmp = tmp(1:2:end)+1i*tmp(2:2:end);
%             INF_Ref2(:,iPE,iCh)=tmp;
%             clear tmp
%         end
%     end
%     %%
%     fclose(fid);
%     
    
    
    vars.ASL.raw.EPI = EPI;
    vars.ASL.raw.refEPI = EPI_Ref;
    vars.ASL.raw.M0 = M0;
    
    vars.PC.raw.KH = PC_KH;
    vars.PC.raw.Ref1 = SUP_Ref1;
    vars.PC.raw.Ref2 = SUP_Ref2;
    
    vars.InfGRE.raw.KH = HbO2_KH;
    vars.InfGRE.raw.Ref1 = INF_Ref1;
    vars.InfGRE.raw.Ref2 = INF_Ref2;
