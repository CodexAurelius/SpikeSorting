

%detect max values
 %determine number of spikes
 %determine threshold value
max_second_col = max(second_col);
descending_second_col = sort(second_col, 'descend');
plot(descending_second_col(1:10));
%disp('top second column values are: '); disp(descending_second_col(1:5));
%disp('max voltage of first electrode is: '); disp(max_second_col);