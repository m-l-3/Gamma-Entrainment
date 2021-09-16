close all; clear all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Initial Parameter Setting %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fs = 250; % Sampling frequency
total_dur = 340; % Duration of entire task in seconds
stim_dur = 40; % Duration of stimulus trials in seconds
rest_dur = 20; % Duration of rest trials in seconds
N = Fs*total_dur; % Total number of samples
trial_num = 6; % Number of stimulus trials

main_data_path = "../Auditory data/";
multi_data_path = "../Multisensory data/";
rest_data_path = "../Rest data/";
part_data_path = "../Partitioned data/";
data_file_info = dir(fullfile(main_data_path, '*.txt'));

Pa_num = length(data_file_info); % Total number of participants
Pa_code = "S"+[1:Pa_num]; % Participant code
% Ignorance list representing whom to exclude from analysis
% S6 and S13 are excluded.
Ignore = repelem("n",Pa_num);
Ignore(6) = "y"; Ignore(13) = "y";
% List of electrodes (10/20 system)
ch_list = ["Fp1","Fp2","F7","F3","Fz","F4","F8","T7","C3","Cz","C4","T8","P7","P3","Pz","P4","P8","O1","O2"];
ch_cell = {"FP1","FP2","F7","F3","FZ","F4","F8","T7","C3","CZ","C4","T8","P7","P3","PZ","P4","P8","O1","O2"};
ch_num = length(ch_list); % Number of electrodes

Fz_idx = find(ch_list == "Fz"); % Index of the Fz channel
Pz_idx = find(ch_list == "Pz"); % Index of the Pz channel

theta_range = [4 8]; % in Hz
gamma_range = [39 41]; % in Hz
adj_40 = [38 42]; % in Hz

addpath('../Code/cbrewer') 
addpath('../Code/plot_topography') 
addpath('../Code/Violinplot-Matlab-master') 

set(groot, ...
'DefaultFigureColor', 'w', ...
'DefaultAxesLineWidth', 1.2, ...
'DefaultAxesXColor', 'k', ...
'DefaultAxesYColor', 'k', ...
'DefaultAxesFontUnits', 'points', ...
'DefaultAxesFontSize', 6, ...
'DefaultAxesFontName', 'Helvetica', ...
'DefaultLineLineWidth', 2, ...
'DefaultTextFontUnits', 'Points', ...
'DefaultTextFontSize', 6, ...
'DefaultTextFontName', 'Helvetica', ...
'DefaultAxesBox', 'off', ...
'DefaultAxesTickLength', [0.02 0.025]);
 
% set the tickdirs to go out - need this specific order
set(groot, 'DefaultAxesTickDir', 'In');
set(groot, 'DefaultAxesTickDirMode', 'manual');

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Occurrence of Entrainment %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grouping the subjects to ent. and non-ent. group

W = 20; % Window length in seconds
nw = Fs*W; % Total number of samples in a window length of W sec.
f = Fs*(0:(nw/2))/nw;
z_th = 3; % Threshold of the z-score values

stim_idx = [1:2; 4:5; 7:8; 10:11; 13:14; 16:17]; % Index of stimulus windows
rest_idx = [3; 6; 9; 12; 15]; % Index of rest windows

freq_idx = f>=adj_40(1) & f<=adj_40(2); % Index of frequecy range aroud 40Hz

load(fullfile(part_data_path, sprintf("f_ch_w%d_Pa.mat", W)))

temp = fData_allWin_allPa(freq_idx, :, stim_idx(:,1), :);
tempz = zscore(abs(temp),0,1);
idx_40 = round(length(tempz)/2); % Index of 40Hz frequency
z_40 = tempz(idx_40,:,:,:);

ent_tr = squeeze(z_40(:, Fz_idx, :, :)) > z_th; % The trial is marked as entrained.
ent = sum(ent_tr) >= trial_num/2+1; % The channel is marked as entrained.
Ent = repelem("n",Pa_num);
Ent(ent) = "y";
Ent

%%
%%%%%% Figs. 1a-b

group = [1:Pa_num] .* (Ignore=="n" & Ent=="n");
group(group==0) = [];

figure
for i=1:length(group)
    subplot(2, 4, i)
    z = squeeze(z_40(:, :, :, group(i)));
    xvalues = {'1','2','3','4','5','6'};
    ent_ch = sum(z > z_th, 2) >= trial_num/2+1; % The channel is marked as entrained.
    ch_copy = ch_list;
    ch_copy(ent_ch==1) = '*  '+ch_copy(ent_ch==1);
    yvalues = ch_copy;
    h = heatmap(xvalues,yvalues,z,'CellLabelColor', 'None', 'fontsize', 5);
    set(struct(h).NodeChildren(3), 'XTickLabelRotation', 45);
    caxis([2 3])
    h.XLabel = Pa_code(group(i));
end

% sg = sgtitle('a','fontweight','bold', 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize', 12);

%%
%%%%%% Fig. 1c

gr_ent = (Ignore=="n" & Ent=="y"); % Index of participants in the ent. group
gr_nonent = (Ignore=="n" & Ent=="n"); % Index of participants in the non-ent. group

z_avg = mean(z_40,3); % Trial average
z_ent = squeeze(mean(z_avg(:, :, :, gr_ent),4)); % Group average
z_nonent = squeeze(mean(z_avg(:, :, :, gr_nonent),4)); % Group average

figure;
subplot(1,2,1)
plot_topography(ch_cell, z_ent)

title('Ent. group','FontSize',6)
colors = cbrewer('seq', 'Blues',9);
colormap(colors);
caxis([min([z_ent,z_nonent]) max([z_ent,z_nonent])])
c = colorbar();c.FontSize=5;

subplot(1,2,2)
plot_topography(ch_cell, z_nonent)

title('Non-ent. group','FontSize',6)
colors = cbrewer('seq', 'Blues',9);
colormap(colors);
caxis([min([z_ent,z_nonent]) max([z_ent,z_nonent])])
c = colorbar();c.FontSize=5;

% sg = sgtitle('c','fontweight','bold', 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize', 12);

%% 
%%%%%% Supp. Fig. 1a
figure
sample_tr = 1; % Sample trail
sample_pa = 1; % Sample participant
temp = fData_allWin_allPa(freq_idx, Fz_idx, sample_tr, sample_pa);
temp = abs(temp);

yyaxis left
plot(f(freq_idx),temp)
y1 = 0; y2 = 0.045;
ylim([y1 y2])

yyaxis right
plot(f(freq_idx),zscore(temp,0,1))
yyaxis right
plot(f(freq_idx),zeros(1,length(f(freq_idx)))+3)
ylim([(y1-mean(temp))/sqrt(var(temp)) (y2-mean(temp))/sqrt(var(temp))])

xlabel('Frequency (Hz)');
yyaxis left
ylabel('Amplitude')
yyaxis right
ylabel('Z-score')
xlim([37.8 42.2])
% sg = sgtitle('a','fontweight','bold', 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize', 12);

%%
%%%%%% Supp. Fig. 1b

x = squeeze(z_avg(:, Fz_idx, :, gr_ent));
y = squeeze(z_avg(:, Fz_idx, :, gr_nonent));

vec = squeeze(z_avg(:, Fz_idx, :, Ignore=="n"));
category = cellstr(Ent(Ignore=="n"));

% Plot violin plots
figure
colors = cbrewer('seq', 'Blues', 9);
v = violinplot(vec, category,'ShowMean',true, 'GroupOrder', {'y','n'},...
    'ViolinColor', colors(8, :), 'ViolinAlpha', 1);
hold on
% Plot error
h = ploterr(1:2,[mean(x) mean(y)], [], [std(x)/sqrt(length(x)) std(y)/sqrt(length(y))], 'k.', 'abshhxy', 0);

% Test of significance
[~, pval, ci, stats] = ttest2(x, y, 'Vartype', 'equal')
mysigstar(gca, [1 2], 8, pval);

set(gca, 'xtick', [1 1.5 2], 'xticklabel', {'Ent.', ch_list(Fz_idx), 'Non-ent.',},'xlim', [0.5 2.5],'XTickLabelRotation',45);
ylabel('Gamma entrainment score'); xlabel('Status');
ylim([-0.5 8.5])
% title('b','FontSize',12, 'FontWeight', 'bold');

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Theta/ Gamma Power analysis %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%% Power time-series (Figs. 2a-b)
range = gamma_range; % or theta-range
W = 40; % window length in seconds
load(fullfile(part_data_path, sprintf("t_ch_w%d_Pa.mat", W)))
stim = squeeze(mean(tData_allWin_allPa(:, Fz_idx, :, :), 3)); % Trial average

W = 20; % window length in seconds
load(fullfile(part_data_path, sprintf("t_ch_w%d_Pa.mat", W)))
rest_idx = [3; 6; 9; 12; 15]; % Index of rest windows
rest = squeeze(mean(tData_allWin_allPa(:, Fz_idx, rest_idx, :), 3)); % Trial average

power_rest = bandpower(rest, Fs, range); % Total rest power

W = 1; % Windows of 1sec
nw = W*Fs; 
w_num = (stim_dur-1)/0.5 + 1; % Number of windows with 0.5sec overlap
powerStim_allWin_allPa = [];
for w=1:w_num
    stim_win = stim((w-1)*nw/2+1:(w-1)*nw/2+nw, :);
    power_stim_win = bandpower(stim_win, Fs, range);
    powerStim_allWin_allPa = cat(1,powerStim_allWin_allPa, power_stim_win./power_rest);
end


Pow_ent = powerStim_allWin_allPa(:, gr_ent);  
Pow_nonEnt = powerStim_allWin_allPa(:, gr_nonent);

figure
colors = cbrewer('seq', 'PuRd', 8);
stdshade(Pow_ent', 0.3, colors(6,:), 1:0.5:stim_dur, 0)
hold on
stdshade(Pow_nonEnt', 0.3, colors(4,:), 1:0.5:stim_dur, 0)

legend('Ent.', 'Non-ent.')
xlabel('Time (s)');
ylabel('Normalized gamma power during stimulation');
xlim([-2 43]);
% ylim([0 4]);
% title('a','FontSize',12);

%%
%%%%% Figs. 2c-d,g
figure

W = 20; % window length in seconds
load(fullfile(part_data_path, sprintf("t_ch_w%d_Pa.mat", W)))

stim_idx = [1:2; 4:5; 7:8; 10:11; 13:14; 16:17]; % Index of stimulus windows
rest_idx = [3; 6; 9; 12; 15]; % Index of rest windows

stim = squeeze(mean(tData_allWin_allPa(:, Fz_idx, stim_idx(:,1), :), 3)); % Trial average
rest = squeeze(mean(tData_allWin_allPa(:, Fz_idx, rest_idx(:,1), :), 3)); % Trial average

range = theta_range; % or gamma_range
power_stim = bandpower(stim, Fs, range);
power_rest = bandpower(rest, Fs, range);
power = power_stim./power_rest; % or power_stim or power_rest

x = power(Ignore=='n' & Ent=='y'); 
y = power(Ignore=='n' & Ent=='n');


range = gamma_range; % or theta_range
power_stim = bandpower(stim, Fs, range);
power_rest = bandpower(rest, Fs, range);
power = power_stim./power_rest; % or power_stim or power_rest

xx = power(Ignore=='n' & Ent=='y'); 
yy = power(Ignore=='n' & Ent=='n');


starloc = 6.5; % Location of stars (the result of test of significance)
colors = cbrewer('seq', 'PuRd', 8);

% Plot violin plots
vec = [x y xx yy];
category = cellstr([repelem("Ent.", length(x)), repelem("Non-ent.", length(y)),...
                    repelem("Ent2.", length(x)), repelem("Non-ent2.", length(y))]);
v = violinplot(vec, category,'ShowMean',true,...
    'GroupOrder', {'Ent.','Non-ent.','Ent2.','Non-ent2.'},...
    'ViolinColor', colors(6, :), 'ViolinAlpha', 1);
hold on
% Plot error
h = ploterr(1:2,[mean(x) mean(y)], [], [std(x)/sqrt(length(x)) std(y)/sqrt(length(y))], 'k.', 'abshhxy', 0);
h = ploterr(3:4,[mean(xx) mean(yy)], [], [std(xx)/sqrt(length(xx)) std(yy)/sqrt(length(yy))], 'k.', 'abshhxy', 0);

% Test of significance
[~, pval, ci, stats] = ttest2(x, y, 'Vartype', 'equal')
mysigstar(gca, [1 2], starloc, pval);
[~, pval, ci, stats] = ttest2(xx, yy, 'Vartype', 'equal')
mysigstar(gca, [3 4], starloc, pval);


set(gca, 'xtick', [1 1.5 2 3 3.5 4], 'xticklabel', {'Ent.', 'Stimulus', 'Non-ent.', 'Ent.', 'Rest', 'Non-ent.'},'xlim', [0.5 4.5],'XTickLabelRotation',45);
% set(gca, 'xtick', [1 1.5 2 3 3.5 4], 'xticklabel', {'Ent.', '\theta', 'Non-ent.', 'Ent.', '\gamma', 'Non-ent.'},'xlim', [0.5 4.5],'XTickLabelRotation',45);
ylabel('Normalized band power'); 
xlabel('Status');

% set(gca,'FontSize',6)
% title('c','FontSize',12)


% [~, pval, ci, stats] = ttest2(x, xx, 'Vartype', 'equal')

%% 
%%%%% Power of theta in 60 sec of Rest data before the main task
% Participants number 9 (S9) and 13 (S13) are missing.
Missing = ["n", "n", "n", "n", "n", "n", "n", "n", "y", "n", "n", "n", "y"];
range = theta_range;

Pa_included = [1:Pa_num] .* (Ignore=="n" & Missing=="n"); % Included participants
Pa_included(Pa_included==0) = [];

% Calculation
power_all = [];
for i=1:length(Pa_included)
  
    tData = load(fullfile(rest_data_path, sprintf('sub-%03d_task-Rest_eeg.txt',Pa_included(i)))); % Data in time domain
    tData_norm = zscore(tData, 0, 1);

    pre_rest = tData_norm(:, Fz_idx);
    
    power = bandpower(pre_rest, Fs, range);
    power_all = cat(1,power_all,power);
end

x = power_all(gr_ent(Pa_included)); 
y = power_all(gr_nonent(Pa_included)); 

%%%%% Fig. 2e
figure()
colors = cbrewer('seq', 'PuRd', 8);
starloc = 0.48;
% Plot violin plots
vec = [x; y];
category = cellstr([repelem("Ent.", length(x)), repelem("Non-ent.", length(y))]);
v = violinplot(vec, category,'ShowMean',true,...
    'GroupOrder', {'Ent.','Non-ent.'},...
    'ViolinColor', colors(6, :), 'ViolinAlpha', 1);
hold on
% Plot error
h = ploterr(1:2,[mean(x) mean(y)], [], [std(x)/sqrt(length(x)) std(y)/sqrt(length(y))], 'k.', 'abshhxy', 0);

% Test of significance
[~, pval, ci, stats] = ttest2(x, y, 'Vartype', 'equal')
mysigstar(gca, [1 2], starloc, pval);

ylabel('Pre-test theta power'); 
xlabel('Status');

% set(gca,'FontSize',6)
% title('e','FontSize',12)
view([90 90])


%% 
%%%%% entanglement of theta to entrainment (Fig. 2f)

range = theta_range;
power_stim = bandpower(stim, Fs, range);

z_avg = mean(z_40,3); % Trial average
% z_avg is named as the gamma entrainment score calculated for each
% participant and each electrode.
z_ent = squeeze(z_avg(:, Fz_idx, :, gr_ent));
z_nonent = squeeze(z_avg(:, Fz_idx, :, gr_nonent));

x1 = z_ent;
x2 = z_nonent;
y1 = power_stim(gr_ent);
y2 = power_stim(gr_nonent);
createFits(x1, y1, x2, y2)

text(mean(x1), max(y1)-mean(y1)/3, sprintf('corr. coeff. = %.2f',corr(x1,y1')))
text(min(x2)+mean(x2)*1.5, max(y2)-mean(y1)/3, sprintf('corr. coeff. = %.2f',corr(x2,y2')))
% title('f','FontSize',12)
% set(gcf, ,'FontSize', 6)

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Temporal Synchrony %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W = 1; % window length in seconds
nw = Fs*W; % Total number of samples in a window length of W sec.
f = Fs*(0:(nw/2))/nw;
load(fullfile(part_data_path, sprintf("f_ch_w%d_Pa.mat", W)))

stim_idx = [1:40 61:100 121:160 181:220 241:280 301:340];
rest_idx = [41:60 101:120 161:180 221:240 281:300];

%%%%% Figs. 3a-b

group = [1:Pa_num] .* (Ignore=="n" & Ent=="y");
group(group==0) = [];

bin = 10; % Number of histogram bins
sample_pa = 1; % Sample participant
ch = [find(ch_list == "Fz") find(ch_list == "Pz")]; % Fz and Pz channels
sz = 36; % Marker size
arrow = 1; % Whether to plot the arrow or not
f40_idx = find(f == 40); % Index of 40Hz frequency component
col = 1; % color. 1 if ent. and 0 if non-ent.

stim_samples = squeeze(fData_allWin_allPa(f40_idx, ch, stim_idx, group(sample_pa)));
rest_samples = squeeze(fData_allWin_allPa(f40_idx, ch, rest_idx, group(sample_pa)));

figure
subplot(2,1,1) % Fz channel

thetap = angle(stim_samples(1,:));
thetaq = angle(rest_samples(1,:));
my_polarhist(thetap, thetaq, bin, arrow, sz, '', col)  
title(ch_list(ch(1)),'FontWeight','normal','FontSize',8);

subplot(2,1,2) % Pz channel

thetap = angle(stim_samples(2,:));
thetaq = angle(rest_samples(2,:));
my_polarhist(thetap, thetaq, bin, arrow, sz, '', col)
title(ch_list(ch(2)),'FontWeight','normal','FontSize',8);

% legend
lgd = legend({'Stimulus', 'Rest'},'Location','east', 'FontSize',8);
hL1 = subplot(2,1,1);
hL2 = subplot(2,1,2);
poshL1 = get(hL1,'position'); 
poshL2 = get(hL2,'position'); 
set(lgd,'position',[poshL1(1)+poshL1(3)*0.5, (poshL1(2)+poshL2(2))*0.48 + poshL2(4)/2, 0.02 0.02]);
title(lgd, 'Ent.')

% sg = sgtitle('a','fontweight','bold', 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize', 12);

%%
%%%%% Figs. 3a appendices

group = [1:Pa_num] .* (Ignore=="n" & Ent=="y");
group(group==0) = [];

bin = 10; % Number of histogram bins
sample_pa = 1; % Sample participant
ch = find(ch_list == "P3"); % Channel
sz = 10; % Marker size
arrow = 1; % Whether to plot the arrow or not
f40_idx = find(f == 40); % Index of 40Hz frequency component
col = 1; % color. 1 if ent. and 0 if non-ent.

stim_samples = squeeze(fData_allWin_allPa(f40_idx, ch, stim_idx, group(sample_pa)));
rest_samples = squeeze(fData_allWin_allPa(f40_idx, ch, rest_idx, group(sample_pa)));

figure

thetap = angle(stim_samples);
thetaq = angle(rest_samples);

my_polarhist(thetap, thetaq, bin, arrow, sz, '', col)  
title(ch_list(ch),'FontWeight','normal','FontSize',8);

%%
% Supp. Figs. 2a-b 

group = [1:Pa_num] .* (Ignore=="n" & Ent=="y");
group(group==0) = [];

bin = 10; % Number of histogram bins
ch = [4 5 6 14 15 16 18 19]; % F3, Fz, F4, P3, Pz, P4, O1 and O2 channels
% ch = [5 15]; % Fz and Pz channels
arrow = 1; % Whether to plot the arrow or not
sz = 5; % Marker size
label = Pa_code(group);
col = 1; % color. 1 if ent. and 0 if non-ent.

figure 
cnt = 0;
for i=1:length(group)
    
    stim_samples = squeeze(fData_allWin_allPa(f40_idx, ch, stim_idx, group(i)));
    rest_samples = squeeze(fData_allWin_allPa(f40_idx, ch, rest_idx, group(i)));

    for j=1:length(ch)

        cnt = cnt+1;
        subplot(length(ch),length(group),(j-1)*length(group)+i) % for ent.
%         subplot(4,4,(j-1)*8+i) % for non-ent.

        thetap = angle(stim_samples(j,:));
        thetaq = angle(rest_samples(j,:));

        if(i==1)
            my_polarhist(thetap, thetaq, bin, arrow, sz, ch_list(ch(j)), col)
        else
            my_polarhist(thetap, thetaq, bin, arrow, sz, '', col)
        end
        
        % for ent.
        if(j==1)
            title(label(i),'FontWeight','normal','FontSize',8);
        end
        % for non-ent.
%         title(label(i),'FontWeight','normal','FontSize',8);
        
        
    end
end

% sg = sgtitle('a','fontweight','bold', 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize', 12);

% % legend  for ent.
lgd = legend({'Stimulus', 'Rest'}, 'FontSize',6);
hL1 = subplot(length(ch),length(group),30);
hL2 = subplot(length(ch),length(group),31);
poshL1 = get(hL1,'position'); 
poshL2 = get(hL2,'position'); 
set(lgd,'position',[poshL1(1) + poshL1(3)*0.85, poshL1(2)-poshL1(4) , 0.1, 0.1]);
title(lgd,'Ent.')

% % legend for Non-ent.
% lgd = legend({'Stimulus', 'Rest'}, 'FontSize',6);
% hL = subplot(4,4,15);
% poshL = get(hL,'position'); 
% set(lgd,'position',[poshL(1)+poshL(3)*1.5, poshL(2)+poshL(4)*0.1 , 0.1, 0.1]);
% title(lgd,'Non-ent.')

%% 
%%%%% Figs. 3c-d

bin = 10; % Number of histogram bins

stim_samples = squeeze(fData_allWin_allPa(f40_idx, :, stim_idx, :));
rest_samples = squeeze(fData_allWin_allPa(f40_idx, :, rest_idx, :));

f = figure('visible','off');
KL = zeros(length(ch), Pa_num); % KL divergence

for ch=1:ch_num
    for i=1:Pa_num
        thetap = angle(stim_samples(ch,:,i));
        thetaq = angle(rest_samples(ch,:,i));

        hp = polarhistogram(thetap,bin);
        pp = hp.Values/length(stim_idx); % Probability distribution of stimulus samples 

        hq = polarhistogram(thetaq,bin);
        qq = hq.Values/length(rest_idx); % Probability distribution of rest samples

        p = pp(pp~=0);
        q = qq(pp~=0);
        kl = -sum(p.*log(q./p)); % KL distance: stim2rest
%         kl = log(bin)+sum(p.*log(p)); % KL distance: stim2unif
        KL(ch,i) = kl;
    end
end

gr_ent = (Ignore=="n" & Ent=="y"); % Index of participants in the ent. group
gr_nonent = (Ignore=="n" & Ent=="n"); % Index of participants in the non-ent. group

x = KL(Fz_idx, gr_ent); % KL distance calculated for Fz channel for the ent. group
xx = KL(Pz_idx, gr_ent); % KL distance calculated for Pz channel for the ent. group
y = KL(Fz_idx, gr_nonent); % KL distance calculated for Fz channel for the non-ent. group
yy = KL(Pz_idx, gr_nonent); % KL distance calculated for Pz channel for the non-ent. group
 
figure()
colors = cbrewer('div', 'BrBG',11);
starloc = 0.85; % Location of stars (the result of test of significance)

% Plot violin plots
vec = [x y xx yy];
category = cellstr([repelem("Ent.", length(x)), repelem("Non-ent.", length(y)),...
                    repelem("Ent2.", length(x)), repelem("Non-ent2.", length(y))]);
v = violinplot(vec, category,'ShowMean',true,...
    'GroupOrder', {'Ent.','Non-ent.','Ent2.','Non-ent2.'},...
    'ViolinColor', colors(10, :), 'ViolinAlpha', 1);
hold on
% Plot error
h = ploterr(1:2,[mean(x) mean(y)], [], [std(x)/sqrt(length(x)) std(y)/sqrt(length(y))], 'k.', 'abshhxy', 0);
h = ploterr(3:4,[mean(xx) mean(yy)], [], [std(xx)/sqrt(length(xx)) std(yy)/sqrt(length(yy))], 'k.', 'abshhxy', 0);

% Test of significance
[~, pval, ci, stats] = ttest2(x, y,'Vartype','equal')
mysigstar(gca, [1 2], starloc, pval);
[~, pval, ci, stats] = ttest2(xx, yy,'Vartype','equal')
mysigstar(gca, [3 4], starloc, pval);

set(gca, 'xtick', [1 1.5 2 3 3.5 4], 'xticklabel', {'Ent.', ch_list(Fz_idx), 'Non-ent.', 'Ent.', ch_list(Pz_idx), 'Non-ent.'},'xlim', [0.5 4.5],'XTickLabelRotation',45);
% ylabel('Stimulus to uniform KL distance');
ylabel('Stimulus to rest KL distance');
xlabel('Status');
% set(gca,'FontSize',6)
% title('c','FontSize',12)

%%
%%%%% Supp. Fig. 2c 
figure

Pval=[];H=[];
for ch=1:ch_num
    x = KL(ch,gr_ent);
    y = KL(ch,gr_nonent);
    [h, pval] = ttest2(x, y,'Vartype','equal');
    H = cat(1,H,h); % h represents whether the result is significant or not
    Pval = cat(1,Pval,pval); % pval represents the p-value
end
new_H = H;
new_H(Pval<0.001) = H(Pval<0.001)+2;
new_H(Pval<0.01 & Pval>0.001) = H(Pval<0.01 & Pval>0.001)+1;

plot_topography(ch_cell, new_H)

colors = cbrewer('div', 'BrBG',11);
colormap(colors(7:10,:));
colorbar('FontSize',6, 'Ticks', 3/8:3/4:3-3/8, 'YTickLabel', ["n.s", "p < 0.05", "p < 0.01", "p < 0.001"])

% sg = sgtitle('c','fontweight','bold', 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize', 12);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Spatial Synchrony %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W = 20; % window length in seconds
load(fullfile(part_data_path, sprintf("t_ch_w%d_Pa.mat", W)))

stim_idx = [1:2; 4:5; 7:8; 10:11; 13:14; 16:17]; % Index of stimulus windows
rest_idx = [3; 6; 9; 12; 15]; % Index of rest windows

stim = squeeze(mean(tData_allWin_allPa(:, :, stim_idx(:,1), :), 3)); % Trial average
rest = squeeze(mean(tData_allWin_allPa(:, :, rest_idx(:,1), :), 3)); % Trial average

H = zeros(ch_num);
Pval = zeros(ch_num);

PLV_diff = zeros(ch_num);

% Calculation of the PLV for each two pair of channels for all participants in
% the gamma range. The statistical measurement would be the difference of PLV
% calculated in the rest and stimulus cycles.

% Running the following commands may take some minutes.
for ch1=1:ch_num % First channel
    for ch2=1:ch_num % Second channel
        if(ch1 == ch2)
            continue
        end
        stim_ch1 = squeeze(stim(:, ch1, :));
        stim_ch2 = squeeze(stim(:, ch2, :));
        rest_ch1 = squeeze(rest(:, ch1, :));
        rest_ch2 = squeeze(rest(:, ch2, :));

        band1 = gamma_range;
        band2 = gamma_range;
        stim_band1_ch1 = angle(hilbert(bandpass(stim_ch1, band1, Fs)));
        stim_band2_ch2 = angle(hilbert(bandpass(stim_ch2, band2, Fs)));
        rest_band1_ch1 = angle(hilbert(bandpass(rest_ch1, band1, Fs)));
        rest_band2_ch2 = angle(hilbert(bandpass(rest_ch2, band2, Fs)));

        z_stim = exp(1i.*(stim_band1_ch1 - stim_band2_ch2));
        z_rest = exp(1i.*(rest_band1_ch1 - rest_band2_ch2));

        PLV_stim = abs(mean(z_stim));
        PLV_rest = abs(mean(z_rest));

        % Statistical measurement
        PLV = PLV_stim - PLV_rest;
        x = PLV(gr_ent);
        y = PLV(gr_nonent);
        
        if(ch1 == Fz_idx && ch2 == Pz_idx) % For Fig. 4a
            x_target = x; y_target = y;
        end
        
        % Test of significance
        [h, pval, ci, stats] = ttest2(x, y,'Vartype','equal');
        H(ch1,ch2) = h;
        Pval(ch1,ch2) = pval;
        
        PLV_diff(ch1,ch2) = mean(x) - mean(y);
    end
end

% We use PLV_diff as the adjacency matrix for plotting the connectivity
% graph (Fig. 4b). We only plot the significant edges.
adj_mat = PLV_diff.*H;
% The adj_mat.edge file would be used in BrainNet Viewer representing the
% edges of the connectivity graph.
save('adj_mat.edge','adj_mat','-ascii');

%%
%%%%% Fig. 4a

figure()
colors = cbrewer('qual', 'Paired', 10);

x = x_target; y = y_target;

% Plot violin plots
vec = [x y];
category = cellstr([repelem("Ent.", length(x)), repelem("Non-ent.", length(y))]);
v = violinplot(vec, category,'ShowMean',true,...
    'GroupOrder', {'Ent.','Non-ent.'},...
    'ViolinColor', colors(6, :), 'ViolinAlpha', 1);
hold on
% Plot error
h = ploterr(1:2,[mean(x) mean(y)], [], [std(x)/sqrt(length(x)) std(y)/sqrt(length(y))], 'k.', 'abshhxy', 0);

mysigstar(gca, [1 2], 0.55, Pval(Fz_idx, Pz_idx));


set(gca, 'xtick', [1 2], 'xticklabel', {'Ent.', 'Non-ent.',},'xlim', [0.5 2.5],'XTickLabelRotation',45);
ylabel('PLV^{Fz Pz}_{stimulus} - PLV^{Fz Pz}_{rest}'); xlabel('Status');
% set(gca, 'FontSize', 6)
% title('a','FontSize',12)
% ylim()


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Theta-Gamma Coupling %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W = 20; % window length in seconds
load(fullfile(part_data_path, sprintf("t_ch_w%d_Pa.mat", W)))

stim_idx = [1:2; 4:5; 7:8; 10:11; 13:14; 16:17]; % Index of stimulus windows
rest_idx = [3; 6; 9; 12; 15]; % Index of rest windows

stim = squeeze(mean(tData_allWin_allPa(:, :, stim_idx(:,1), :), 3)); % Trial average
rest = squeeze(mean(tData_allWin_allPa(:, :, rest_idx(:,1), :), 3)); % Trial average

%%%%% Calculating MVL comodulograms (for two joint single frequencies in 
% theta/gamma range) for included participants for the Fz channel.  

% Running the following commands may take several minutes. You can skip the
% following commands in this part since the MVL is calculated before and
% saved here.

Pa_included = [1:Pa_num] .* (Ignore=="n"); % Included participants
Pa_included(Pa_included==0) = [];

highfreq = 35:45; % High frequency range for amplitude
lowfreq = 4:12; % Low frequency range for phase
amp_length = length(highfreq);
phase_length = length(lowfreq);

% Time-frequency MVL
tf_MVL_all = zeros(amp_length, phase_length, 2, Pa_num);

% Running the following loops may take several minutes
tic;
for sub=1:length(Pa_included)
    sub % subject

    for win=[1 2] % 1 for stimulus and 2 for rest windows     
        if(win==1)
            x = stim(:, Fz_idx, Pa_included(sub));
        else
            x = rest(:, Fz_idx, Pa_included(sub));
        end
        for i = 1:phase_length
            for j = 1:amp_length
              l_freq = lowfreq(i);
              h_freq = highfreq(j);
              tf_MVL_all(j, i, win, Pa_included(sub)) = tfMVL(x, h_freq, l_freq, Fs);
            end
        end
    end
end
toc;
save("Amp_Phase_StimRest_allPa_Fz.mat",'tf_MVL_all');

%% 
%%%%% Fig. 5a
figure()

highfreq = 35:45; % High frequency range for amplitude
lowfreq = 4:12; % Low frequency range for phase

PAC = load('Amp_Phase_StimRest_allPa_Fz.mat');
tf_MVL_all = PAC.tf_MVL_all;
mi = min(min(min(min(tf_MVL_all(:, :, :, [1 2])))));
ma = max(max(max(max(tf_MVL_all(:, :, :, [1 2])))));

ylab = ["Stimulus", "rest"];
tit = ["Ent.", "Non-ent."];
colors = cbrewer('seq', 'Greens',20);

for i=1:2
    for j=1:2
        subplot(2,2,(i-1)*2+j)
            plot_comodulogram(tf_MVL_all(:, :, i, j), highfreq, lowfreq)
            caxis([mi ma])
            if(j==1)
                ylabel(ylab(i),'FontSize',6)
            end
            xticks([4 6 8 10 12])
            if(i==1)
                title(tit(j),'FontSize',6);
            end
    end
end

colormap(colors)

hL = subplot(2,2,4);
pos = get(hL,'position');
colorbar('Position', [pos(1)+pos(3)+0.02  pos(2)+pos(4)*0.25 0.01  pos(4)*1.5]);
text(2.6,34,'Amplitude frequncy (Hz)','Rotation',90, 'FontSize',6);
text(0.85,33,'Phase frequncy (Hz)','FontSize',6);

% sg = sgtitle('a','fontweight','bold', 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize', 12);

%%
%%%%% Supp. Figs. 3a-b

group = [1:Pa_num] .* (Ignore=="n" & Ent=="y");
group(group==0) = [];

PAC = load('Amp_Phase_StimRest_allPa_Fz.mat');
tf_MVL_all = PAC.tf_MVL_all;

Pa_included = [1:Pa_num] .* (Ignore=="n");
Pa_included(Pa_included==0) = [];

mi = min(min(min(min(tf_MVL_all(:, :, :, Pa_included)))));
ma = max(max(max(max(tf_MVL_all(:, :, :, Pa_included)))));

figure
colors = cbrewer('seq', 'Greens',50);

for sub=1:length(group)

    subplot(length(group), 2, (sub-1)*2+1)
    
    plot_comodulogram(tf_MVL_all(:, :, 1, group(sub)), highfreq, lowfreq) % Stimulus
    
    xticks([4 6 8 10 12])
    caxis([mi*2 ma/2])
    ylabel(Pa_code(group(sub)),'FontSize',6)
    if(sub==length(group))
        xlabel('Stimulus','FontSize',6)
    end
    
    hL = subplot(length(group), 2, sub*2);
    pos = get(hL,'position');
    
    plot_comodulogram(tf_MVL_all(:, :, 2, group(sub)), highfreq, lowfreq) % Rest
    
    xticks([4 6 8 10 12])
    caxis([mi*2 ma/2])
    if(sub==length(group))
        xlabel('Rest','FontSize',6)
    end
    
end

text(2.7,32,'Amplitude frequncy (Hz)','Rotation',90);
text(1.2,30,'Phase frequncy (Hz)');

colormap(colors);
colorbar('Position', [pos(1)+pos(3)+0.02  pos(2)+pos(4) 0.01  pos(4)*3]); % for ent. group
% colorbar('Position', [pos(1)+pos(3)+0.02  pos(2)+pos(4)*1.6 0.01  pos(4)*6]); % for non-ent. group

% sg = sgtitle('a','fontweight','bold', 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize', 12);

%%
%%%%% Calculating MVL in frequency bands for all participants and all channels

% Running the following commands may take some minutes. You can skip this
% part since the MVL is calculated before and saved here.

high = gamma_range; % High frequency range for amplitude
low = theta_range; % Low frequency range for phase
highfreq = high(1):1:high(2);
amp_length = length(highfreq);
lowfreq = low(1):1:low(2);
phase_length = length(lowfreq);

mvl_band = zeros(ch_num, 2, Pa_num);

% Running the following loops may take some minutes.
tic;
for sub=13:13
    for win=1:2 % 1 for stimulus and 2 for rest windows
        for ch=1:ch_num
            if (win==1)
                x = stim(:, ch, sub);
            else
                x = rest(:, ch, sub);
            end  
            mvl_band(ch, win, sub) = band_tfMVL(x, highfreq, lowfreq, Fs);
        end
    end
end
toc;

save("Ch_StimRest_allPa_bandMVL.mat",'mvl_band');

%%
%%%%% Figs. 5b, d

PAC_band = load('Ch_StimRest_allPa_bandMVL.mat');
mvl_band = PAC_band.mvl_band;

mvl_band_stim = squeeze(mvl_band(:, 1, :));
mvl_band_rest = squeeze(mvl_band(:, 2, :));
mvl_band_diff = mvl_band_stim - mvl_band_rest;

% mvl_measure = mvl_band_diff; % for Fig. 5b
mvl_measure = mvl_band_stim; % for Fig. 5d

x = mvl_measure(Fz_idx, gr_ent);
y = mvl_measure(Fz_idx, gr_nonent);

xx = mvl_measure(Pz_idx, gr_ent);
yy = mvl_measure(Pz_idx, gr_nonent);

starloc = 3.3; % Location of stars (the result of test of significance)
colors = cbrewer('qual', 'Paired', 8);

figure
% Plot violin plots
vec = [x y xx yy];
category = cellstr([repelem("Ent.", length(x)), repelem("Non-ent.", length(y)),...
                    repelem("Ent2.", length(x)), repelem("Non-ent2.", length(y))]);
v = violinplot(vec, category,'ShowMean',true,...
    'GroupOrder', {'Ent.','Non-ent.','Ent2.','Non-ent2.'},...
    'ViolinColor', colors(4, :), 'ViolinAlpha', 1);
hold on
% Plot error
h = ploterr(1:2,[mean(x) mean(y)], [], [std(x)/sqrt(length(x)) std(y)/sqrt(length(y))], 'k.', 'abshhxy', 0);
h = ploterr(3:4,[mean(xx) mean(yy)], [], [std(xx)/sqrt(length(xx)) std(yy)/sqrt(length(yy))], 'k.', 'abshhxy', 0);

% Test of significance
[~, pval, ci, stats] = ttest2(x, y,'Vartype','equal')
mysigstar(gca, [1 2], starloc, pval);
[~, pval, ci, stats] = ttest2(xx, yy,'Vartype','equal')
mysigstar(gca, [3 4], starloc, pval);


set(gca, 'xtick', [1 1.5 2 3 3.5 4], 'xticklabel', {'Ent.', 'Fz', 'Non-ent.', 'Ent.', 'Pz', 'Non-ent.'},'xlim', [0.5 4.5],'XTickLabelRotation',45,'FontSize',6);
xlabel('Status');
% ylabel('MVL_{stimulus} - MVL_{rest}'); % for Fig. 5b
ylabel('MVL_{stimulus}'); % for Fig. 5d
% set(gca, 'FontSize',6)
% title('d', 'FontSize',12)


%%
%%%%% Fig. 5c
figure()

coup1 = mean(mvl_band_diff(:, gr_ent), 2); % Group average
coup2 = mean(mvl_band_diff(:, gr_nonent), 2); % Group average

coup12 = cat(2,coup1,coup2); 
mi = min(min(coup12));
ma = max(max(coup12));
colors = cbrewer('seq', 'Greens',8);

subplot(1,2,1)
plot_topography(ch_cell, coup1)

title('Ent. group','FontSize',6)
caxis([mi ma])

subplot(1,2,2)
plot_topography(ch_cell, coup2)

title('Non-ent. group','FontSize',6)
caxis([mi ma])

colormap(colors);
colorbar('FontSize',6)

% sg = sgtitle('c','fontweight','bold', 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize', 12);

%%
%%%%% Supp. Fig. 3c
figure()

H=[];
Pval=[];

for ch=1:ch_num
    x = mvl_band_diff(ch, gr_ent);
    y = mvl_band_diff(ch, gr_nonent);
    [h, pval] = ttest2(x, y,'Vartype','equal');
    H = cat(1,H,h); % h represents whether the result is significant or not
    Pval = cat(1,Pval,pval); % pval represents the p-value
end

new_H = H;
new_H(Pval<0.001) = H(Pval<0.001)+2;
new_H(Pval<0.01 & Pval>0.001) = H(Pval<0.01 & Pval>0.001)+1;

plot_topography(ch_cell, new_H)

colors = cbrewer('seq', 'Greens',4);
colormap(colors);
colorbar('FontSize',6)
colorbar('FontSize',6, 'Ticks', 3/8:3/4:3-3/8, 'YTickLabel', ["n.s", "p < 0.05", "p < 0.01", "p < 0.001"])

% sg = sgtitle('c','fontweight','bold', 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize', 12);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Audio + Visual Entrainment %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% New parameter setting
Fs = 500; % Sampling frequency
total_dur = 40; % Duration of each epoch in seconds
N = Fs*total_dur; % Total number of samples
epoch_num = 3; % Number of epochs

multi_data = load(fullfile(multi_data_path, "sub-1_task-MultisensoryGammaEntrainment_eeg.txt"));
% This data consists of three epochs of 40sec EEG recorded signal during
% 1)auditory 2)visual and 3)audio-visual stimulations respectively.

auditory = multi_data(1:N, :);
visual = multi_data(N+1:2*N, :);
audio_visual = multi_data(2*N+1:end, :);

epoched_data = cat(3, auditory, visual, audio_visual);
epoched_data_norm = zscore(epoched_data,0,1); % Normalized data

% Using the first 20sec of each epoch
main_data = epoched_data_norm(1:20*Fs, :, :);

%%
%%%%% Calculating MVL comodulograms (for two joint single frequencies in 
% theta/gamma range) for the Fz channel.

% Running the following commands may take several minutes. You can skip the
% following, since the MVL is calculated before and saved.

highfreq = 35:45; % High frequency range for amplitude
lowfreq = 4:8; % Low frequency range for phase
amp_length = length(highfreq);
phase_length = length(lowfreq);

% Time-frequency MVL
tf_MVL_all = zeros(amp_length, phase_length, epoch_num);
% Running the following loops may take several minutes
for epoch=1:epoch_num
    x = main_data(:, Fz_idx, epoch);
    tic;
    for i = 1:phase_length
        for j = 1:amp_length
          l_freq = lowfreq(i);
          h_freq = highfreq(j);
          [tf_MVL_all(j, i, epoch)] = tfMVL(x, h_freq, l_freq, Fs);
        end
    end
    toc;
end

save("Amp_Phase_AudVisAudvis_Fz.mat",'tf_MVL_all');

%%
%%%%%% Fig. 6

PAC = load('Amp_Phase_AudVisAudvis_Fz.mat');
tf_MVL_all = PAC.tf_MVL_all;

highfreq = 35:45; % High frequency range for amplitude
lowfreq = 4:8; % Low frequency range for phase
amp_length = length(highfreq);
phase_length = length(lowfreq);

ma = max(max(max(tf_MVL_all)))/2;
mi = min(min(min(tf_MVL_all)));
colors = cbrewer('seq', 'Oranges',50);

figure
subplot(1,3,1);

plot_comodulogram(tf_MVL_all(:,:,1), highfreq, lowfreq)
caxis([mi ma])
xlabel('Phase frequency (Hz)', 'Fontsize', 6);ylabel('Amplitude frequency (Hz)', 'Fontsize', 6)
title('Auditory', 'Fontsize', 6)
set(gca,'FontSize',6);

subplot(1,3,2);

plot_comodulogram(tf_MVL_all(:,:,2),highfreq,lowfreq)
caxis([mi ma])
xlabel('Phase frequency (Hz)', 'Fontsize', 6);ylabel('Amplitude frequency (Hz)', 'Fontsize', 6)
title('Visual', 'Fontsize', 6)
set(gca,'FontSize',6);

subplot(1,3,3);

plot_comodulogram(tf_MVL_all(:,:,3),highfreq,lowfreq)
caxis([mi ma])
xlabel('Phase frequency (Hz)', 'Fontsize', 6);ylabel('Amplitude frequency (Hz)', 'Fontsize', 6)
title('Auditory + Visual', 'Fontsize', 6)
set(gca,'FontSize',6);

colormap(colors);
hL = subplot(1,3,3);
pos = get(hL,'position');
colorbar('Position', [pos(1)+pos(3)+0.05  pos(2)+pos(4)/24 0.01  0.75]);


%%%%%%%%%%%%%%%%%%%%%
%%%%%% The End %%%%%%
%%%%%%%%%%%%%%%%%%%%%
