
clc
clear
root_dir="C:\Users\sinad\OneDrive - Emory University\Patient Temp\R1638E_RepFR1_Session1\R1638E_2022-10-17_17-35-39";
% openning the file 
file_path = fullfile(root_dir , "eeg_data_reduced_chA.edf");
if isfile(file_path)
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
L = size(temp_data);
t = (0:L(1)-1)*T(1);
figure(1)
subplot(2,2,1)
plot(t,temp_data)
xlabel("time (sec)")
ylabel("Amp (uV)")

% params.Fs = data_info.samples/ data_info.duration;
params.fpass = [0.1 500];
params.err = [1 0.05];
params.trialave = 0;
params.tapers=[3 5];
% Plotting PSD using mtspectrumc funciton
[S,f,Serr]=mtspectrumc(temp_data, params); 
subplot(2,2,2)
plot(f,log10(S));
xlabel("freq (Hz)")
ylabel("Power (log)")

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
subplot(2,2,4)
plot(f_filt,log10(S_filt));
xlabel("freq (Hz)")
ylabel("Power (log)")

%% Preprocess the iEEG recording and a template: 2015 Horak method
f_down = round(params.Fs/200);
X_raw = downsample(temp_data_filt,f_down);
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
subplot(2,2,2)
plot([1,2,3,4],kernel)
xlabel("x")
ylabel("y")

subplot(2,2,3)
plot(t_preproc, X_preproc)
xlabel("time (sec)")
ylabel("Amp (uV)")
title("conv")

subplot(2,2,4)
plot(t_preproc(1:end-1), diff(X_preproc))
xlabel("time (sec)")
ylabel("Amp (uV)")
title("diff")
%% 2009 Nonclercq preprocessing method: 
% instead of convolution - derivative > squaring > moving window integration

X_preproc_diff = diff(X_raw);
X_preproc_sqr = sign(X_preproc_diff).*power(X_preproc_diff,2);
X_preproc_smooth = movmean(X_preproc_sqr,10);
figure(3)
subplot(4,1,1)
plot(t_preproc, X_raw)
xlabel("time (sec)")
ylabel("Amp (uV)")
title("raw")

subplot(4,1,2)
plot(t_preproc(1:end-1), X_preproc_diff)
xlabel("time (sec)")
ylabel("Amp (uV)")
title("diff")

subplot(4,1,3)
plot(t_preproc(1:end-1), X_preproc_sqr)
xlabel("time (sec)")
ylabel("Amp (uV)")
title("sqr")

subplot(4,1,4)
plot(t_preproc(1:end-1), X_preproc_smooth)
xlabel("time (sec)")
ylabel("Amp (uV)")
title("movmean")