close all; clear all; clc

Fs = 250; % Sampling frequency
total_dur = 340; % Duration of entire main task in seconds
stim_dur = 40; % Duration of stimulus trials in seconds
rest_dur = 20; % Duration of rest trials in seconds
N = Fs*total_dur; % Total number of samples

main_data_path = "../Auditory data/";
part_data_path = "../Partitioned data/";
data_file_info = dir(fullfile(main_data_path, '*.txt'));

Pa_num = length(data_file_info); % Total number of participants
Pa_code = "S"+[1:Pa_num]; % Participant code

% List of electrodes (10/20 system)
ch_list = {"Fp1","Fp2","F7","F3","Fz","F4","F8","T7","C3","Cz","C4","T8","P7","P3","Pz","P4","P8","O1","O2"};
ch_num = length(ch_list);

%% Auditory Data partitioning 

% Set W to desired window length
W = 20; % in seconds
% Number of windows
w_num = total_dur/W; % or 6 for W=40. 
nw = Fs*W; % Total number of samples in a window length of W sec.
f = Fs*(0:(nw/2))/nw;

tData_allWin_allPa = [];
fData_allWin_allPa = [];

for i=1:Pa_num
    tData = load(fullfile(main_data_path, sprintf('sub-%03d_task-AuditoryGammaEntrainment_eeg.txt',i))); % Data in time domain
    tData = tData(1:N,:);
    tData_norm = zscore(tData,0,1); % Normalizing the data
    tData_allWin = [];
    fData_allWin = [];
    
    for w=1:w_num
        tData_win = tData_norm( (w-1)*nw+1:w*nw,: ); % tData in window length of W sec.
%         tData_win = tData_norm( (w-1)*nw*3/2+1: (w-1)*nw*3/2 +nw,: );% For 40sec trials of stimulus only

        % Calculating FFT
        fData_win = fft(tData_win);
        P2 = (fData_win/nw);
        P1 = P2(1:nw/2+1,:);
        P1(2:end-1,:) = 2*P1(2:end-1,:);
        
        fData_allWin = cat(3,fData_allWin,P1);
        tData_allWin = cat(3,tData_allWin,tData_win);
    end
    tData_allWin_allPa = cat(4,tData_allWin_allPa,tData_allWin);
    fData_allWin_allPa = cat(4,fData_allWin_allPa,fData_allWin);
end

save(fullfile(part_data_path, sprintf("t_ch_w%d_Pa.mat", W)),'tData_allWin_allPa');
save(fullfile(part_data_path, sprintf("f_ch_w%d_Pa.mat", W)),'fData_allWin_allPa');

%%%%% Epoched data are 4-D arrays containing:
% 1st dimention: Time/Frequency samples of each window
% 2nd dimention: Channels
% 3rd dimention: Windows 
% 4th dimention: Participants
