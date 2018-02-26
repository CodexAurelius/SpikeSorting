clear all;
close all;

% Input file

% This code loads the "EMG_example_2_fs_2k.csv" file and assigns the
% individual channel we would like to test. You have to extend the code so
% that the other two files "EMG_example_1_90s_fs_2k.csv" and
% "EMG_example_3_fs_10k" can be loaded for doing the channel assignment. Spike sorting is
% done only for the selected channel. The channel will be randomly selected
% for the evaluation so that you have to make sure your algorithm works
% well for all channels.


    %M= csvread('EMG_example_1_90s_fs_2k.csv'); %read in csv file
%M= csvread('EMG_example_3_fs_10k.csv'); %read in csv file

M= csvread('EMG_example_2_fs_2k.csv'); %read in csv file
time= M(:,1); % first column is the time series
fs= (time(2)-time(1))^-1; % calculate the sample frequecy
channel_number= size(M,2)-1; % num of channels in the database
for i=1:channel_number,
figure('Color',[1 1 1]);plot(time,M(:,i+1)); %plot each channel
str= sprintf('Channel %d',i);
xlabel('seconds');title(str);xlim([time(1) time(size(time,1))]); % label and title each plots
end
channel_select= 1; % select channel for testing. channel_select <= channel_number
test_input= M(:,channel_select+1); % only this test_input will go through the spike sorting algorithm
