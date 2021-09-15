% Author: Chang Hi Lee
% Title: Applying Perceptron Classification Model to Classify Locust Olfactory Data
% Description: This code attempts to classify locust olfactory response to two types of odors, Hexanol(pink) and Benzaldehyde(blue). 
% The Perceptron classifier(a.k.a "filter" in code) is trained on two sets of responses from pink and blue stimuli, respectively. 
% Then, it predicts the response to a new set of pink and blue stimuli. The accuracy comparisons of each prediction are plotted. 

clear all
close all
load Data.mat
%Analysis 1
figure
hold on
plot(a1_time-1.5, a1_voltageresponse+1.9,'linewidth',3)
plot(a1_time-1.5, a1_stim./50,'linewidth',3)
xlim([0 22])
title('Hexanol EAG Response to 1s Stimulus with 2s ISI')
xlabel('Time (s)'); ylabel('Voltage (mV)')
legend('EAG Response (mV)', 'Valve Open/Close (Unitless)')
set(gca,'fontsize',30)

%Analysis 2
a1_vfilt = filtfilt(ones(1,200)/200,1,a1_voltageresponse);
% 200 point smoothing (assuming 1000 samples/s sampling rate)
figure
hold on
plot(a1_time-1.5, a1_vfilt+2.2,'linewidth',3)
plot(a1_time-1.5, a1_stim./50,'linewidth',3)
xlim([0 22])
title('Hexanol EAG Response to 1s Stimulus with 2s ISI, 200 Point Smoothing')
xlabel('Time (s)'); ylabel('Voltage (mV)')
legend('EAG Response (mV)', 'Valve Open/Close (Unitless)')
set(gca,'fontsize',30)

figure
a2_2_vfilt = filtfilt(ones(1,200)/200,1,a2_2_vresponse);
plot(a2_2_time-1.7, a2_2_vfilt+2.5, a2_2_time-1.7, a2_2_stim./50,'linewidth',3)
xlim([0 20]); ylim([-0.4 1.2])
title('Hexanol EAG Response to 1s Stimulus with 1s ISI, 200 Point Smoothing')
xlabel('Time (s)'); ylabel('Voltage (mV)')
legend('EAG Response (mV)', 'Valve Open/Close (Unitless)')
set(gca,'fontsize',30)

figure
a2_3_vfilt = filtfilt(ones(1,200)/200,1,a2_3_vresponse);
plot(a2_3_time-1, a2_3_vfilt+2.45, a2_3_time-1, a2_3_stim./50,'linewidth',3)
xlim([0 15.2])
title('Hexanol EAG Response to 1s Stimulus with 500ms ISI, 200 Point Smoothing')
xlabel('Time (s)'); ylabel('Voltage (mV)')
legend('EAG Response (mV)', 'Valve Open/Close (Unitless)')
set(gca,'fontsize',30)

figure
a2_4_vfilt = filtfilt(ones(1,200)/200,1,a2_4_vresponse);
plot(a2_4_time-2.3, a2_4_vfilt+2.1, a2_4_time-2.3, a2_4_stim./50,'linewidth',3)
xlim([0 14.5])
title('Hexanol EAG Response to 1s Stimulus with 250ms ISI, 200 Point Smoothing')
xlabel('Time (s)'); ylabel('Voltage (mV)')
legend('EAG Response (mV)', 'Valve Open/Close (Unitless)')
set(gca,'fontsize',30)

figure
a2_5_vfilt = filtfilt(ones(1,200)/200,1,a2_5_vresponse);
plot(a2_5_time-1.8, a2_5_vfilt+2, a2_5_time-1.8, a2_5_stim./50,'linewidth',3)
% xlim([0 14.5])
title('Hexanol EAG Response to Random Stimulus Length and ISI, 200 Point Smoothing')
xlabel('Time (s)'); ylabel('Voltage (mV)')
legend('EAG Response (mV)', 'Valve Open/Close (Unitless)')
set(gca,'fontsize',30)

figure
bluea2_1_vfilt = filtfilt(ones(1,200)/200,1,bluea2_1_vresponse);
plot(bluea2_1_time-3, bluea2_1_vfilt+2, bluea2_1_time-3, bluea2_1_stim./50,'linewidth',3)
xlim([0 29])
title('Benzaldehyde EAG Response to 1s Stimulus with 2s ISI, 200 Point Smoothing')
xlabel('Time (s)'); ylabel('Voltage (mV)')
legend('EAG Response (mV)', 'Valve Open/Close (Unitless)')
set(gca,'fontsize',30)

figure
bluea2_2_vfilt = filtfilt(ones(1,200)/200,1,bluea2_2_vresponse);
plot(bluea2_2_time-1, bluea2_2_vfilt+2.425, bluea2_2_time-1, bluea2_2_stim./50,'linewidth',3)
xlim([0 20])
title('Benzaldehyde EAG Response to 1s Stimulus with 1s ISI, 200 Point Smoothing')
xlabel('Time (s)'); ylabel('Voltage (mV)')
legend('EAG Response (mV)', 'Valve Open/Close (Unitless)')
set(gca,'fontsize',30)

figure
bluea2_3_vfilt = filtfilt(ones(1,200)/200,1,bluea2_3_vresponse);
plot(bluea2_3_time-1.1, bluea2_3_vfilt+2.83, bluea2_3_time-1.1, bluea2_3_stim./50,'linewidth',3)
xlim([0 17]); ylim([-0.2 1])
title('Benzaldehyde EAG Response to 1s Stimulus with 500ms ISI, 200 Point Smoothing')
xlabel('Time (s)'); ylabel('Voltage (mV)')
legend('EAG Response (mV)', 'Valve Open/Close (Unitless)')
set(gca,'fontsize',30)

figure
bluea2_4_vfilt = filtfilt(ones(1,200)/200,1,bluea2_4_vresponse);
plot(bluea2_4_time-1, bluea2_4_vfilt+2.66, bluea2_4_time-1, bluea2_4_stim./50,'linewidth',3)
xlim([0 12.7]); ylim([-0.2 1])
title('Benzaldehyde EAG Response to 1s Stimulus with 250ms ISI, 200 Point Smoothing')
xlabel('Time (s)'); ylabel('Voltage (mV)')
legend('EAG Response (mV)', 'Valve Open/Close (Unitless)')
set(gca,'fontsize',30)

figure
bluea2_5_vfilt = filtfilt(ones(1,200)/200,1,bluea2_5_vresponse);
plot(bluea2_5_time-1, bluea2_5_vfilt+2.68, bluea2_5_time-1, bluea2_5_stim./50,'linewidth',3)
xlim([0 16.5]); ylim([-0.2 0.8])
title('Benzaldehyde EAG Response to Random Stimulus Length and ISI, 200 Point Smoothing')
xlabel('Time (s)'); ylabel('Voltage (mV)')
legend('EAG Response (mV)', 'Valve Open/Close (Unitless)')
set(gca,'fontsize',30)

%pink
Stim_row15 = Stimulus(15,:);
Response_row15 = Response(15,:);
down_Stim_row15 = downsample(Stimulus(15,:),10);
down_Response_row15 = downsample(Response(15,:),10);

%blue
Stim_row16 = Stimulus(16,:);
Response_row16 = Response(16,:);
down_Stim_row16 = downsample(Stimulus(16,:),10);
down_Response_row16 = downsample(Response(16,:),10);


S_pink = zeros(5000, 500);
S_pink(1,:) = down_Stim_row15(1:500);
j=1;
for i=2:5000
    S_pink(i,:) = down_Stim_row15(j+1:j+500);
    if j<4500
    j=j+1;
    else
    end
end
S_pink = [S_pink ones(5000,1)];

W_pink = pinv(S_pink)*down_Response_row15';

% Stim_pink = downsample(a2_5_stim(381:15380),3);
% Stim_blue = downsample(bluea2_5_stim(565:15564),3);
% S_pink = zeros(5000, 500);
S_blue = zeros(5000, 500);
S_blue(1,:) = down_Stim_row16(1:500);
%for blue
m=1;
for i=2:5000
    S_blue(i,:) = down_Stim_row16(m+1:m+500);
    if m<4500
    m=m+1;
    else
    end
end
S_blue = [S_blue ones(5000,1)];
W_blue = pinv(S_blue)*down_Response_row16';

R_pink = S_pink*W_pink;
R_blue = S_blue*W_blue;

time=0.005:0.005:25;

R_pink_filt = filtfilt(ones(1,200)/200,1,R_pink);
R_blue_filt = filtfilt(ones(1,200)/200,1,R_blue);

figure
hold on
plot(time, filtfilt(ones(1,200)/200,1,down_Response_row15)-4,'linewidth',3)
plot(time, R_pink_filt-4,'linewidth',3)
title('Hexanol Actual and Model Stimulus Response')
xlabel('Time (s)'); ylabel('Voltage (mV)')
legend('Actual Response','Model Response')
set(gca,'fontsize', 30)

figure
hold on
plot(time, filtfilt(ones(1,200)/200,1,down_Response_row16)-9,'linewidth',3)
plot(time,R_blue_filt-9,'linewidth',3)
title('Benzaldehyde Actual and Model Stimulus Response')
xlabel('Time (s)'); ylabel('Voltage (mV)')
legend('Actual Response','Model Response')
set(gca,'fontsize', 30)

subplot(1,2,1)
plot(-5:0.01:0,W_pink,'linewidth',3)
xlim([-6 1]); ylim([-2 12])
title('Hexanol Filter')
xlabel('Time (s)'); ylabel('Gain (unitless)')
set(gca,'fontsize',30)
subplot(1,2,2)
plot(-5:0.01:0,W_blue,'linewidth',3)
xlim([-6 1])
title('Benzaldehyde Filter')
xlabel('Time (s)'); ylabel('Gain (unitless)')
set(gca,'fontsize',30)


% create S(stimulus) matrix from Row 15 of sample data
% R is the response, W is the transfer function, we want to get the transfer
% function(i.e. "train" the data) that accurately captures the behavior of
% the locusts(can also use the two random ISI data sets that we did as trainer data)
% After obtaining the transfer function, apply it to the experimental data
% and compare the results to see if the transfer function is accurate.
% y = downsample([Row 15], 10) (because Row 15 is 1000 points long)


% ISI : 2s, 1s, 500ms, 250ms -> Freq : 
% pink : [0.995, 0.93, 0.799, 0.716]  (mV)
% blue : [0.769, 0.603, 0.594, 0.589]  (mV)
freq = 1./[2, 1, 0.5, 0.25];
amp_pink = [0.995, 0.93, 0.799, 0.716];
amp_blue = [0.769, 0.603, 0.594, 0.589];
f_pink=fit(freq',amp_pink','smoothingspline')
f_blue=fit(freq',amp_blue','smoothingspline')

figure
hold on
plot(f_pink, freq, amp_pink)
plot(f_blue, freq, amp_blue)
xlim([0 4.5])
title('Frequency Response of Olfactory System')
xlabel('Frequency (Hz)'); ylabel('Voltage (mV)')
set(gca,'fontsize',30)















