
clc
clear
root_dir="C:\Users\sinad\OneDrive - Emory University\Patient Temp\R1638E_RepFR1_Session1\R1638E_2022-10-17_17-35-39";
% openning the file 
file_path = fullfile(root_dir , "eeg_data_reduced_chA.edf");
if isfile(file_path)
    addpath("C:\Users\sinad\OneDrive - Georgia Institute of Technology\DrGross\Lou\Analysis_CommonCode\downloadedFunctions\edfread_v2.10.0.1\");
    [data_info, data] = edfread(file_path); 
% data_info = edfinfo(file_path);
end

% nCh = data_info.NumSignals;
time_seg = 10; %for a 10 sec segment
params.Fs = data_info.samples(1)/data_info.duration;

Num_samples=time_seg*params.Fs;
temp_data = zeros(round(Num_samples),1);
for i=1:Num_samples+1
        temp_data(i)=data(1,i+round(10*params.Fs));
end
% temp_data={};
% for j = 1: nCh
% 
%     % Plotting 10-20 sec of data
%     
%     for i=1:101
%         temp_data{i,j}=[data{i,j}{:}];
%     end
% end
% github example:
% https://github.com/Willie-Lab/Analysis_FlickerProject/blob/bb191838d685c2822714c2d83fac4c216ffb1cec/Analysis/run_PSD_flicker.m 

T = 1/params.Fs;
[L,~] = size(temp_data);
t = (0:L-1)*T(1);
figure(1)
subplot(2,2,1)
plot(t,temp_data)
xlabel("time (sec)")
ylabel("Amp (uV)")
title("x raw")

% params.Fs = data_info.samples/ data_info.duration;
params.fpass = [0.1 500];
params.err = [1 0.05];
params.trialave = 0;
params.tapers=[3 5];
% Plotting PSD using mtspectrumc funciton
[S,f,~]=mtspectrumc(temp_data, params); 

Y = fft(temp_data);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f2 = params.Fs*(0:(L/2))/L;

subplot(2,2,2)
plot(f2,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1 (f)|')


% plot(f,log10(S));
% xlabel("freq (Hz)")
% ylabel("Power (log)")
% title("X freq spectrogram")

% Filter line noise

d = designfilt('bandstopiir','FilterOrder',4, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',params.Fs);

temp_data_filt = filtfilt(d, temp_data);
[S_filt,f_filt,Serr]=mtspectrumc(temp_data_filt, params);

subplot(2,2,3)
plot(t,temp_data_filt)
xlabel("time (sec)")
ylabel("Amp (uV)")
title("x filterd 60 Hz line noise")
subplot(2,2,4)
plot(f_filt,log10(S_filt));
xlabel("freq (Hz)")
ylabel("Power (log)")
title("X freq spectrogram")
%% IED detection: 2015 Horak method
% Step 1: Preprocess the iEEG recording and a template 

f_down = round(params.Fs/200);
X_raw = downsample(temp_data_filt,f_down,0);
T_preproc = 1/200;
L_preproc = size(X_raw);
t_preproc = (0:L_preproc(1)-1)*T_preproc(1);

kernel = [-2,-1,1,2]/8;
X_preproc = conv(X_raw,kernel,'same');


figure(2)
subplot(2,2,1)
plot(t_preproc, X_raw)
xlabel("time (sec)")
ylabel("Amp (uV)")
title("X raw downsampled to 200 Hz")

subplot(2,2,2)
plot([1,2,3,4],kernel)
xlabel("x")
ylabel("y")
title("kernel")

subplot(2,2,3)
plot(t_preproc, X_preproc)
xlabel("time (sec)")
ylabel("Amp (uV)")
title("conv(Xraw, kernel)")

subplot(2,2,4)
plot(t_preproc(1:end-1), diff(X_raw))
xlabel("time (sec)")
ylabel("Amp (uV)")
title("diff(Xraw)")
%% 2009 Nonclercq preprocessing method: 
% instead of convolution - derivative > squaring > moving window integration

X_preproc_diff = diff(X_raw);
X_preproc_sqr = sign(X_preproc_diff).*power(X_preproc_diff,2);
X_preproc_smooth = movmean(X_preproc_sqr,10);
% figure(3)
% subplot(4,1,1)
% plot(t_preproc, X_raw)
% xlabel("time (sec)")
% ylabel("Amp (uV)")
% title("raw")
% 
% subplot(4,1,2)
% plot(t_preproc(1:end-1), X_preproc_diff)
% xlabel("time (sec)")
% ylabel("Amp (uV)")
% title("diff")
% 
% subplot(4,1,3)
% plot(t_preproc(1:end-1), X_preproc_sqr)
% xlabel("time (sec)")
% ylabel("Amp (uV)")
% title("sqr")
% 
% subplot(4,1,4)
% plot(t_preproc(1:end-1), X_preproc_smooth)
% xlabel("time (sec)")
% ylabel("Amp (uV)")
% title("movmean")
%% Step2: covolve the signal with a triangular template
triangle_win = triang(200*.06);
template = conv(triangle_win,kernel,"same");

X_ccorr = conv(X_preproc, template, "same" );

figure(3)
subplot(2,2,1)
plot(1:12, triangle_win)
xlabel("sample")
xlim([1 12])
ylabel("y")
title("triangular window")

subplot(2,2,2)
plot(1:12, template)
xlabel("sample")
xlim([1 12])
ylabel("y")
title("template = conv(triangleWin, kernel)")

subplot(2,2,3)
plot(t_preproc, X_preproc)
xlabel("time (sec)")
ylabel("Amp (uV)")
title("Xconv = conv(Xraw, kernel)")

subplot(2,2,4)
plot(t_preproc, X_ccorr)
xlabel("time (sec)")
ylabel("Amp (uV)")
title("Xccorr = conv(Xconv, template)")
%% Step 3: Normalize the cross-corrolated by the median std dev. from 1 sec sliding window
% 1 sec sliding window: for 200 Hz signal is a 200 sample wide window
mov_std_dev = movstd(X_ccorr,200);
med_mov_stdDev = median(mov_std_dev);
X_ccorr_norm = normalize(X_ccorr,'scale',med_mov_stdDev);

%% Step 4: label local peaks above an emperical threshold
signal_findpeaks_path = "C:\Program Files\MATLAB\R2022a\toolbox\signal\signal\";
if exist(signal_findpeaks_path,"dir")
    addpath(signal_findpeaks_path);
else
    addpath("C:\Program Files\MATLAB\R2021b\toolbox\signal\signal")
end
std_dev_norm = std(X_ccorr_norm);
[X_pks, X_loc] = findpeaks(X_ccorr_norm,"MinPeakHeight",3*std_dev_norm); % , "MinPeakHeight",2.5

figure(4)
subplot(2,2,1)
plot(t_preproc, mov_std_dev)
xlabel("time (sec)")
ylabel("Amp (uV)")
title("moving std dev")

subplot(2,2,2)
plot(t_preproc, X_ccorr_norm)
xlabel("time (sec)")
ylabel("Amp (uV)")
title("Xccorr norm = normalize(Xccorr, median stdDev)")

subplot(2,2,3)
hold on
plot(t_preproc, X_ccorr_norm)
scatter(X_loc/200,X_pks)        
xlabel("time (sec)")
ylabel("Amp (uV)")
title("Xccorr norm with peaks annotated")
hold off