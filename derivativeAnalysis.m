

%derivative analysis
der_second_col = M(numel(second_col));
for idx= 2:numel(second_col)
    der_second_col(idx) = second_col(idx)-second_col(idx-1); 
end
clc;
plot(first_col, der_second_col);
axis([0,21,-1,1]);

%find average of plot
%find threshold of "spike"
%where do zeros and spikes correspond to on unmodified signal