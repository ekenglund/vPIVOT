%% obtain quantitative information from perfusion time course - from image
close all,clear all,clc
a = dir('20241103*.mat');
load(a.name)
%%
Perf_img = varsASL.Perf.Img;
Perf_tc = varsASL.Perf.AverageTimecourse_IMG_mean;
figure, plot(Perf_tc), 
title('select min and max baseline index')
[x,~] = ginput(2);

for k = 1:7
    Perf_tc_norm(:,k) = Perf_tc(:,k) - mean(Perf_tc(floor(x(1)):floor(x(2)),k));
end
clear x
figure, plot(Perf_tc_norm), pause(1);

%% next step is to define period of exercise
for k = 1:180
    temp_img = Perf_img(:,:,k);
    size_nan(k) = length(find(isnan(temp_img)));
end

%assume that lots of nans are present at the start and end of exercise

figure, plot(size_nan), hold on, plot(Perf_tc_norm(:,1))
title('Define period of exercise: click the start of exercise')
[x(1), y(1)] = ginput(1);
hold on, plot(x(1),y(1),'r.','MarkerSize',20)
plot(x(1)+23,y(1),'r.', 'MarkerSize',20)

title('Define period of exercise: click the end of exercise')
[x(2),~] = ginput(1);
ex_start_time = 4*floor(x(1));
ex_start_idx = floor(x(1));
ex_end_time = 4*ceil(x(2));
ex_end_idx = ceil(x(2));
% start_search_idx = ceil(x(3));
%% perfusion time course is smoothed
Perf_tc_norm_ex_end = zeros(length(ex_end_idx:180),7);
for k = 1:7
Perf_tc_norm_ex_end(:,k) = smooth(Perf_tc_norm(ex_end_idx:end,k));
% Perf_tc_norm_ex_end(:,k) = smooth(Perf_tc_norm_ex_end(:,k));
end
%%
clear x y
figure, for k = 1:7
    plot(Perf_tc_norm_ex_end(:,k)), title('Define region to search for peak (click start and end search area). If you want to exclude, click>100'), [x(:,k),y(:,k)]=ginput(2); end

end_search_formax = floor(min(x(2,:)));
start_search_formax = round(x(1,:));

for k = 1:7
    Perf_tc_norm_ex_searchformax(:,k) = Perf_tc_norm_ex_end(1:end_search_formax,k);
end
%%
Peak_Perfusion = zeros(1,7);
TTP = zeros(1,7);
TTP_idx = zeros(1,7);
AUC_idx = zeros(1,7);
TTR = zeros(1,7);
AUC = zeros(1,7);
AUC_1_2_idx = zeros(1,7);
TTR_1_2 = zeros(1,7);
AUC_1_2 = zeros(1,7);

for k = 1:7
    if start_search_formax(k)<100
        Peak_Perfusion(k) = max(Perf_tc_norm_ex_searchformax(start_search_formax(k):end,k));
        TTP(k) = find(Perf_tc_norm_ex_searchformax(:,k)==Peak_Perfusion(k))*4;
        TTP_idx(k) = TTP(k)/4;
%         AUC_idx(k) = find(Perf_tc_norm_ex_end(:,k)<0,1)-1;
        AUC_idx(k) = find(Perf_tc_norm_ex_end(TTP_idx(k):end,k)<0,1)+TTP_idx(k)-1;
        TTR(k) = AUC_idx(k)*4;
        AUC(k) = sum(Perf_tc_norm_ex_end(start_search_formax:AUC_idx(k),k))*4;
        AUC_1_2_idx(k) = find(Perf_tc_norm_ex_end(TTP_idx(k):end,k)<Peak_Perfusion(k)/2,1)-1+TTP_idx(k);
        TTR_1_2(k) = AUC_1_2_idx(k)*4;
        AUC_1_2(k) = sum(Perf_tc_norm_ex_end(start_search_formax:AUC_1_2_idx(k),k))*4;

        figure, plot(Perf_tc_norm_ex_end(:,k)), hold on, plot(TTP_idx(k),Peak_Perfusion(k),'r.','MarkerSize',20), plot(AUC_idx(k),0,'.','MarkerSize',20), plot(AUC_1_2_idx(k),Peak_Perfusion(k)/2,'.','MarkerSize',20), drawnow, title(k), pause(1); hold off,
    else
        disp(['skipping idx = ',num2str(k)]),
    end
end

Export_params = [Peak_Perfusion,TTP,TTR,AUC, TTR_1_2, AUC_1_2];



%%

save('20241104_vPIVOT3_30_ExportParams.mat')
% close all,clear all,clc
% load('20240104_ExportParams.mat');
% 
% % Peak_Perfusion = zeros(1,7);
% % TTP = zeros(1,7);
% % TTP_idx = zeros(1,7);
% AUC_1_2_idx = zeros(1,7);
% TTR_1_2 = zeros(1,7);
% AUC_1_2 = zeros(1,7);
% for k = 1:7
%     if start_search_formax(k)<100
% %         Peak_Perfusion(k) = max(Perf_tc_norm_ex_searchformax(start_search_formax(k):end,k));
% %         TTP(k) = find(Perf_tc_norm_ex_searchformax(:,k)==Peak_Perfusion(k))*4;
% %         TTP_idx(k) = TTP(k)/4;
%         AUC_1_2_idx(k) = find(Perf_tc_norm_ex_end(TTP_idx(k):end,k)<Peak_Perfusion(k)/2,1)-1+TTP_idx(k);
%         TTR_1_2(k) = AUC_1_2_idx(k)*4;
%         AUC_1_2(k) = sum(Perf_tc_norm_ex_end(start_search_formax:AUC_1_2_idx(k),k))*4;
%         figure, plot(Perf_tc_norm_ex_end(:,k)), hold on, plot(TTP_idx(k),Peak_Perfusion(k),'r.','MarkerSize',20), plot(AUC_idx(k),0,'.','MarkerSize',20), plot(AUC_1_2_idx(k),Peak_Perfusion(k)/2,'.','MarkerSize',20), drawnow, title(k), pause(1); hold off,
%     else
%         disp(['skipping idx = ',num2str(k)]),
%     end
% end
% 
% Export_params = [Peak_Perfusion,TTP,TTR,AUC, TTR_1_2, AUC_1_2];
% 
% save('20240104_ExportParams.mat')
