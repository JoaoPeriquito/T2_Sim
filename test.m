%%%%%%%%T2(s) simulator%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Description%%%%%%%%%%%%%%%%%
% for a giving range of T2(*) and SNR values it calculates the associated  fitting error% 

%%%%%%%%%%%%%

clc
clear all
close all
tic

%% INPUT %%
freq_MR = 128;
%%%%%%%%%%%%%%

out_phase = (1/(freq_MR*3.5).*1000).*[1:40];
first_echo = out_phase(2);
threshold_error = 4; %plots will be displayed with the absolute error in T2s from [0 threshold_error]


%% T2s Parameters
T2s_values    = [10:2:30];
TEs_used      = [first_echo out_phase(2).*[2:40]];
Max_Te = find(TEs_used>=max(T2s_values).*3);
TEs_used = TEs_used(1:Max_Te(1));
M0              = 1;

numOfSimulations = 1;

SI = zeros(length(T2s_values), length(TEs_used));

for i=1:length(T2s_values)
SI(i,:) = M0.*exp(-TEs_used./T2s_values(i));
end

%% Add Noise
SNR = [10 11 12 13 14 16 18 20 23 26 30 35 40 45 50 60 70 80 100 140 200 400 1000];
sigma1 = zeros (length(T2s_values),length(SNR));

for i=1:length(T2s_values)
sigma1(i,:) = SI(i,1)./SNR;
end

ImgTest = zeros(length(T2s_values),length(SNR),numOfSimulations,length(TEs_used));

for k=1:numOfSimulations
    %% Create Syntetic Data
    for i=1:length(SNR)
        for j =1:length(T2s_values)
        ImgTest(j,i,k,:) = abs(squeeze(SI(j,:)) + sigma1(j,i) * randn(size(TEs_used)));
        end
    end
end

for k=1:numOfSimulations
    for i=1:length(SNR)
        for j =1:length(T2s_values)
        [M0_start, T2s_start] = startValues (TEs_used, squeeze(ImgTest(j,i,k,:)));
        [num2str(T2s_values(j)) '_' num2str(T2s_start)];
        
%Longest T2*
        cuttof_longest = find(TEs_used>=max(T2s_values));
        %         TEs_used(cuttof_longest(1))
        [x,ssq,cnt] = fitT2s(TEs_used(1:cuttof_longest(1)), squeeze(ImgTest(j,i,k,1:cuttof_longest(1))),M0_start,T2s_start);
       
        fit_results.longest.M0{j,i,k} =x(1);
        fit_results.longest.t2s{j,i,k} =x(2);
        fit_results.longest.ssq{j,i,k} =ssq;
        fit_results.longest.cnt{j,i,k} =cnt;
        clear x ssq cnt
        
%cuttof T2* real        
        cuttof_real = find(TEs_used>=T2s_values(j));
%         TEs_used(cuttof_real(1))
        [x,ssq,cnt] = fitT2s(TEs_used(1:cuttof_real(1)), squeeze(ImgTest(j,i,k,1:cuttof_real(1))),M0_start,T2s_start);
        fit_results.cuttofReal.M0{j,i,k} =x(1);
        fit_results.cuttofReal.t2s{j,i,k} =x(2);
        fit_results.cuttofReal.ssq{j,i,k} =ssq;
        fit_results.cuttofReal.cnt{j,i,k} =cnt;
        clear x ssq cnt
% 
% %cuttof T2* approx        
        cuttof_approx = find(TEs_used>=T2s_start);
%         TEs_used(cuttof_approx(1))
        [x,ssq,cnt] = fitT2s(TEs_used(1:cuttof_real(1)), squeeze(ImgTest(j,i,k,1:cuttof_real(1))),M0_start,T2s_start);
        fit_results.cuttof_approx.M0{j,i,k} =x(1);
        fit_results.cuttof_approx.t2s{j,i,k} =x(2);
        fit_results.cuttof_approx.ssq{j,i,k} =ssq;
        fit_results.cuttof_approx.cnt{j,i,k} =cnt;
        clear x ssq cnt
%         
% %cuttof T2* real * 1.5
         cuttof_real_1p5 = find(TEs_used>=T2s_values(j).*1.5);
%         TEs_used(cuttof_real_1p5(1))
        [x,ssq,cnt] = fitT2s(TEs_used(1:cuttof_real_1p5(1)), squeeze(ImgTest(j,i,k,1:cuttof_real_1p5(1))),M0_start,T2s_start);
        fit_results.cuttofReal1p5.M0{j,i,k} =x(1);
        fit_results.cuttofReal1p5.t2s{j,i,k} =x(2);
        fit_results.cuttofReal1p5.ssq{j,i,k} =ssq;
        fit_results.cuttofReal1p5.cnt{j,i,k} =cnt;
        clear x ssq cnt

%cuttof T2* approx * 1.5        
        cuttof_approx_1p5 = find(TEs_used>=T2s_start.*1.5);
%         TEs_used(cuttof_approx(1))
        [x,ssq,cnt] = fitT2s(TEs_used(1:cuttof_approx_1p5(1)), squeeze(ImgTest(j,i,k,1:cuttof_approx_1p5(1))),M0_start,T2s_start);
        fit_results.cuttof_approx_1p5.M0{j,i,k} =x(1);
        fit_results.cuttof_approx_1p5.t2s{j,i,k} =x(2);
        fit_results.cuttof_approx_1p5.ssq{j,i,k} =ssq;
        fit_results.cuttof_approx_1p5.cnt{j,i,k} =cnt;
        clear x ssq cnt
        
 %cuttof T2* real * 2.0
         cuttof_real_2p0 = find(TEs_used>=T2s_values(j).*2.0);
%         TEs_used(cuttof_real_2p0(1))
        [x,ssq,cnt] = fitT2s(TEs_used(1:cuttof_real_2p0(1)), squeeze(ImgTest(j,i,k,1:cuttof_real_2p0(1))),M0_start,T2s_start);
        fit_results.cuttofReal2p0.M0{j,i,k} =x(1);
        fit_results.cuttofReal2p0.t2s{j,i,k} =x(2);
        fit_results.cuttofReal2p0.ssq{j,i,k} =ssq;
        fit_results.cuttofReal2p0.cnt{j,i,k} =cnt;
        clear x ssq cnt
        
%cuttof T2* approx * 2.0
        cuttof_approx_2p0 = find(TEs_used>=T2s_start.*2.0);
%         TEs_used(cuttof_approx_2p0(1))
        [x,ssq,cnt] = fitT2s(TEs_used(1:cuttof_approx_2p0(1)), squeeze(ImgTest(j,i,k,1:cuttof_approx_2p0(1))),M0_start,T2s_start);
        fit_results.cuttof_approx_2p0.M0{j,i,k} =x(1);
        fit_results.cuttof_approx_2p0.t2s{j,i,k} =x(2);
        fit_results.cuttof_approx_2p0.ssq{j,i,k} =ssq;
        fit_results.cuttof_approx_2p0.cnt{j,i,k} =cnt;
        clear x ssq cnt
        
% %cuttof T2* real * 2.5
% 
%          cuttof_real_2p5 = find(TEs_used>=T2s_values(j).*2.5);
% %         TEs_used(cuttof_real_2p5(1))
%         [x,ssq,cnt] = fitT2s(TEs_used(1:cuttof_real_2p5(1)), squeeze(ImgTest(j,i,k,1:cuttof_real_2p5(1))),M0_start,T2s_start);
%         fit_results.cuttofReal2p5.M0{j,i,k} =x(1);
%         fit_results.cuttofReal2p5.t2s{j,i,k} =x(2);
%         fit_results.cuttofReal2p5.ssq{j,i,k} =ssq;
%         fit_results.cuttofReal2p5.cnt{j,i,k} =cnt;
%         clear x ssq cnt
% 
% %cuttof T2* approx * 2.5
%         
%         cuttof_approx_2p5 = find(TEs_used>=T2s_start.*2.5);
% %         TEs_used(cuttof_approx(1))
%         [x,ssq,cnt] = fitT2s(TEs_used(1:cuttof_approx_2p5(1)), squeeze(ImgTest(j,i,k,1:cuttof_approx_2p5(1))),M0_start,T2s_start);
%         fit_results.cuttof_approx_2p5.M0{j,i,k} =x(1);
%         fit_results.cuttof_approx_2p5.t2s{j,i,k} =x(2);
%         fit_results.cuttof_approx_2p5.ssq{j,i,k} =ssq;
%         fit_results.cuttof_approx_2p5.cnt{j,i,k} =cnt;
%         clear x ssq cnt
        end
    end
end
toc

%Average the iterations
for i=1:length(SNR)
    for j =1:length(T2s_values)
        error_t2s_longest(j,i) = abs((T2s_values(j) - mean([fit_results.longest.t2s{j,i,:}])));
        error_t2s_cuttofReal(j,i) = abs((T2s_values(j) - mean([fit_results.cuttofReal.t2s{j,i,:}])));
        error_t2s_cuttof_approx(j,i) = abs((T2s_values(j) - mean([fit_results.cuttof_approx.t2s{j,i,:}])));
        error_t2s_cuttofReal1p5(j,i) = abs((T2s_values(j) - mean([fit_results.cuttofReal1p5.t2s{j,i,:}])));
        error_t2s_cuttof_approx_1p5(j,i) = abs((T2s_values(j) - mean([fit_results.cuttof_approx_1p5.t2s{j,i,:}])));
        error_t2s_cuttofReal2p0(j,i) = abs((T2s_values(j) - mean([fit_results.cuttofReal2p0.t2s{j,i,:}])));
        error_t2s_cuttof_approx_2p0(j,i) = abs((T2s_values(j) - mean([fit_results.cuttof_approx_2p0.t2s{j,i,:}])));
%         error_t2s_cuttofReal2p5(j,i) = abs((T2s_values(j) - mean([fit_results.cuttofReal2p5.t2s{j,i,:}])));
%         error_t2s_cuttof_approx_2p5(j,i) = abs((T2s_values(j) - mean([fit_results.cuttof_approx_2p5.t2s{j,i,:}])));


    end
end

figure, imagesc(error_t2s_longest),caxis([0 threshold_error]), axis image,ylabel('true t2*'),xlabel('SNR'), title(['freq: ' num2str(freq_MR) ' error t2s longest']),set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23] ), set(gca,'XTickLabel',SNR ),set(gca,'YTick',[1 2 3 4 5 6 7 8 9 10 11] ), set(gca,'YTickLabel',T2s_values ),hcb=colorbar; title(hcb,'absolute error in T2*')
figure, imagesc(error_t2s_cuttofReal),caxis([0 threshold_error]), axis image,ylabel('true t2*'),xlabel('SNR'), title(['freq: ' num2str(freq_MR) ' error t2s cuttofReal']),set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23] ), set(gca,'XTickLabel',SNR ),set(gca,'YTick',[1 2 3 4 5 6 7 8 9 10 11] ), set(gca,'YTickLabel',T2s_values ),hcb=colorbar; title(hcb,'absolute error in T2*')
figure, imagesc(error_t2s_cuttof_approx),caxis([0 threshold_error]), axis image,ylabel('true t2*'),xlabel('SNR'), title(['freq: ' num2str(freq_MR) ' error t2s cuttof approx']),set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23] ), set(gca,'XTickLabel',SNR ),set(gca,'YTick',[1 2 3 4 5 6 7 8 9 10 11] ), set(gca,'YTickLabel',T2s_values ),hcb=colorbar; title(hcb,'absolute error in T2*')
figure, imagesc(error_t2s_cuttofReal1p5),caxis([0 threshold_error]), axis image,ylabel('true t2*'),xlabel('SNR'), title(['freq: ' num2str(freq_MR) ' error t2s cuttofReal1p5']),set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23] ), set(gca,'XTickLabel',SNR ),set(gca,'YTick',[1 2 3 4 5 6 7 8 9 10 11] ), set(gca,'YTickLabel',T2s_values ),hcb=colorbar; title(hcb,'absolute error in T2*')
figure, imagesc(error_t2s_cuttof_approx_1p5),caxis([0 threshold_error]), axis image,ylabel('true t2*'),xlabel('SNR'), title(['freq: ' num2str(freq_MR) ' error t2s cuttof approx 1p5']),set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23] ), set(gca,'XTickLabel',SNR ),set(gca,'YTick',[1 2 3 4 5 6 7 8 9 10 11] ), set(gca,'YTickLabel',T2s_values ),hcb=colorbar; title(hcb,'absolute error in T2*')
figure, imagesc(error_t2s_cuttofReal2p0),caxis([0 threshold_error]), axis image,ylabel('true t2*'),xlabel('SNR'), title(['freq: ' num2str(freq_MR) ' error t2s cuttofReal2p0']),set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23] ), set(gca,'XTickLabel',SNR ),set(gca,'YTick',[1 2 3 4 5 6 7 8 9 10 11] ), set(gca,'YTickLabel',T2s_values ),hcb=colorbar; title(hcb,'absolute error in T2*')
figure, imagesc(error_t2s_cuttof_approx_2p0),caxis([0 threshold_error]), axis image,ylabel('true t2*'),xlabel('SNR'), title(['freq: ' num2str(freq_MR) ' error t2s cuttof approx 2p0']),set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23] ), set(gca,'XTickLabel',SNR ),set(gca,'YTick',[1 2 3 4 5 6 7 8 9 10 11] ), set(gca,'YTickLabel',T2s_values ),hcb=colorbar; title(hcb,'absolute error in T2*')
% figure, imagesc(error_t2s_cuttofReal2p5),caxis([0 threshold_error]), axis image,ylabel('true t2*'),xlabel('SNR'), title(['freq: ' num2str(freq_MR) ' error t2s cuttofReal2p5']), set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23] ), set(gca,'XTickLabel',SNR ), set(gca,'YTickLabel',T2s_values ),hcb=colorbar; title(hcb,'absolute error in T2*')
% figure, imagesc(error_t2s_cuttof_approx_2p5),caxis([0 threshold_error]), axis image,ylabel('true t2*'),xlabel('SNR'), title(['freq: ' num2str(freq_MR) ' error t2s cuttof approx 2p5']), set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23] ), set(gca,'XTickLabel',SNR ), set(gca,'YTickLabel',T2s_values ),hcb=colorbar; title(hcb,'absolute error in T2*')
% 
