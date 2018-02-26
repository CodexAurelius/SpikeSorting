%Load CSV Data
filename1 = 'EMG_example_3_90s_fs_2k.csv';
filename2 = 'EMG_example_2_fs_2k_.csv';
filename3 = 'EMG_example_3_fs_10k.csv';
M = csvread(filename3);
first_col = M(:,1); %this is time
second_col = M(:,2); %electrode #1
third_col = M(:,3); %electrode #2
%plot(M(:,1));
%plot(M(:,2));
plot(first_col, second_col);
%disp('row 2 column 1 value is: '); disp(M(2,1));





    
    
    




%clear all;
%clc;

%Open data file
%fid = fopen('EMG_example_3_fs_10k.csv');
%readData = textscan(fid,'%s%f%f','delimiter',',');
%fclose(fid);
%date = datevec(readData{1});
%col1 = readData{1};
%col2 = readData{2};
%Read data in CSV file
%readData = textscan(fid,'%f %f %f', 'HeaderLines', 1, 'Delimiter ',',');



%Extract data from read data
%error, perhaps because no .csv file???
%xData = readData(1,1)(:,1);
%y1Data = readData(1,2)(:,2);
%y2Data = readData(1,3)(:,1);

%Plot data
%f1 = figure(1);
%cla; hold on; grid on;
%plot(col1, col2, 'k-');
%plot(col1, col2,'x');
%plot(xData, y2Data, 'r-');

%%Read Header LIne

%Open data file


%Read data from CSV file
%readHeader = textscan(fid, %s
