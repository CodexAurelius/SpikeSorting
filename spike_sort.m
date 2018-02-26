%load CSV data channel by channel
%M= csvread('EMG_example_2_fs_2k.csv'); %read in csv file
M= csvread('EMG_example_3_fs_10k.csv'); %read in csv file
%M= csvread('EMG_example_1_90s_fs_2k.csv'); %read in csv file
    time= M(:,1); % first column is the time series
    fs= (time(2)-time(1))^-1; % calculate the sample frequecy
    channel_number= size(M,2)-1; % num of channels in the database
       for i=1:channel_number,
           figure('Color',[1 1 1]);plot(time,M(:,i+1)); %plot each channel
           str= sprintf('Channel %d',i);
           xlabel('seconds');title(str);xlim([time(1) time(size(time,1))]); % label and title each plots
       end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Spike Data in Tabular Format');
 for j = 1:channel_number
         channel_select = j; % select channel for testing. channel_select <= channel_number
         test_input = M(:,channel_select+1); % only this test_input will go through the spike sorting algorithm
         
 %filter signal via numerical derivative
         der_input = der_filter_signal(test_input, time);
            %plot each channel
            %figure('Color',[1 1 1]);plot(time,der_input); str= sprintf('Derivative %d',j); xlabel('seconds');title(str); xlim([time(1) time(size(time,1))]); 
       
 %detect spikes
         num_spikes_guess = 65;
         df = detect_max_base(der_input,num_spikes_guess);
         rdf = remove_overlap_detect_max_base(der_input,num_spikes_guess);
          %for k = 1:4,   figure('Color',[1 1 1]); plot(der_input(df(k,2): df(k,3)) ); str= sprintf('First Spike %d-%d',j,k);  xlabel('time');title(str);%xlim([time(1) time(size(time,1))]); % label and title each plot
           %end 

 %align spikes
         %center peak voltages at zero
%         figure('Color',[1 1 1]);
%          for k = 1:num_spikes_guess
%               plot(der_input(df(k,1)-10: df(k,1)+10) );
%               %plot(der_input(df(k,1)-(df(k,1) - df(k,2)): df(k,1) + (df(k,3) - df(k,1))) );
%               %plot(der_input(df(k,2): df(k,3) ));
%               hold on
%               str= sprintf('Overlapped spikes');
%                 xlabel('time');title(str);%xlim([time(1) time(size(time,1))]); % label and title each plot
%          end
         
            figure(8+j);
            %disp(size(rdf,1));
          %figure('Color',[1 1 1]);
          for k = 1:size(rdf,1)
              %plot(der_input(df(k,1)-(df(k,1) - df(k,2)): df(k,1) + (df(k,3) - df(k,1))) );
              plot(der_input(rdf(k,1)-20: rdf(k,1)+20) );
                %plot(der_input(df(k,2): df(k,3) ));
              hold on
              str= sprintf('Aligned spikes: Channel %d',j);
               xlabel('time');title(str);%xlim([time(1) time(size(time,1))]); % label and title each plot
          end
         
 %extract features
    %feature extraction occures during remove_overlap_detect_max_base, and more features
    %are extracting before clustering  inextract_features_cluster
 %cluster spikes

        %plots 1x1 clusters
       features = extract_features_cluster(rdf, der_input);
       
 %classify spikes
        fprintf('Index of spike max, Index of start, Index of end, Variance, PCA1, PCA2, Mean voltage');
        disp(rdf);
        if(j<channel_number) %dont print on last iteration
            fprintf('Spike Data in Tabular Format, continued:\n');
        end
        if(j==channel_number)
            fprintf('End of Spike Data in Tabular Format');
        end
        

 end
 
 
 








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function der_filtered_input = der_filter_signal(input, time)
    %derivative analysis
    der_col = zeros(1,numel(input)); %row matrix
    for idx= 2:numel(input)
        der_col(idx) = (input(idx)-input(idx-1))/(time(idx) - time(idx-1)); 
    end
    der_filtered_input = der_col;
%     clc;
%     plot(first_col, der_col);
%     axis([0,21,-1,1]);
end

function filtered_input = filter_signal(input)
    %derivative analysis
    row = zeros(1,numel(input)); %row matrix
    for idx= 2:numel(input)
        row(1,idx) = input(idx);
    end
    filtered_input = row;
%     clc;
%     plot(first_col, der_col);
%     axis([0,21,-1,1]);
end

function [descending_indexes_of_local_maxs_bases] = detect_max_base(input,num_spikes_guess)
    descending_indexes_of_local_maxs = zeros(numel(input),1 ); %col matrix
    descending_indexes_of_local_maxs_bases = zeros(size(input,1),3);
    base_inds = zeros(numel(input),2); %column matrix
    input = abs(input); %no negative spikes
    descending_col = sort(input, 'descend'); %sort magnitude of voltage data, 
    %disp(descending_col(2));

    %do for all local maxs
    for j = 1:num_spikes_guess
        descending_indexes_of_local_maxs(j,1) = find(input == descending_col(j),1); %find index in input of highest value
        loc_max_ind = descending_indexes_of_local_maxs(j,1);
        bases = find(~input); %list of indexes of input zeros
        bases_right = find(bases >= loc_max_ind); %find closest index of zero to the right
        bases_left = find(bases < loc_max_ind); %find closest index of zero to the left
        base_inds(j,2) = min( abs(  loc_max_ind - bases(bases_right(1):bases_right(size(bases_right,2))))) +loc_max_ind;
        base_inds(j,1) = -1* min( abs(  loc_max_ind - bases(bases_left(1):bases_left(size(bases_left,2))))) +loc_max_ind;
        %remove submat of baseinds
        %if element found in descending_col matches elements from base to
        %base of spike, skip
        
        descending_indexes_of_local_maxs_bases(j,1) = descending_indexes_of_local_maxs(j,1); %maxes to first col
        descending_indexes_of_local_maxs_bases(j,2) = base_inds(j,1); %left base to second col
        descending_indexes_of_local_maxs_bases(j,3) = base_inds(j,2); %right base to third col
        
        
%         fprintf('max value is %d\n', descending_col(j)); %debugging
%         fprintf('index of max value is %d\n', loc_max_ind); %debugging
%          fprintf('index of left base value is %d\n', base_inds(j,1)); %debugging
%          fprintf('index of right base value is %d\n\n', base_inds(j,2)); %debugging
          %figure(j+7); plot(input(base_inds(j,1):base_inds(j,2)));
            %fprintf('min value is %d\n', descending_col(size(descending_col))); %debugging
    end
        
end

function [descending_indexes_of_local_maxs_bases] = remove_overlap_detect_max_base(input,num_spikes_guess)
    num_features = 7;
    descending_indexes_of_local_maxs = zeros(numel(input),1 ); %col matrix
    descending_indexes_of_local_maxs_bases = zeros(size(input,1),num_features);
    base_inds = zeros(numel(input),2); %column matrix
    input = abs(input); %no negative spikes
    descending_col = sort(input, 'descend'); %sort magnitude of voltage data, 
    %halfwidth = zeros(numel(input),1 ); %col matrix
    %disp(descending_col(1,1:10));

    %do for all local maxs
    c = 1;
    for j = 1:num_spikes_guess
        if(find(input == descending_col(1,j),1))
            descending_indexes_of_local_maxs(c,1) = find(input == descending_col(1,j),1); %find index in input of highest value
            
            loc_max_ind = descending_indexes_of_local_maxs(c,1);
            bases = find(~input); %list of indexes of input zeros
            bases_right = find(bases >= loc_max_ind); %find closest index of zero to the right
            bases_left = find(bases < loc_max_ind); %find closest index of zero to the left
            base_inds(c,2) = min( abs(  loc_max_ind - bases(bases_right(1):bases_right(size(bases_right,2))))) +loc_max_ind;
            base_inds(c,1) = -1* min( abs(  loc_max_ind - bases(bases_left(1):bases_left(size(bases_left,2))))) +loc_max_ind;
            
            %extract features
            descending_indexes_of_local_maxs_bases(c,1) = descending_indexes_of_local_maxs(c,1); %maxes to first col
            descending_indexes_of_local_maxs_bases(c,2) = base_inds(c,1); %left base to second col
            descending_indexes_of_local_maxs_bases(c,3) = base_inds(c,2); %right base to third col
            descending_indexes_of_local_maxs_bases(c,4) = var(input(base_inds(c,1):base_inds(c,2))); %variance
            descending_indexes_of_local_maxs_bases(c,7) = mean(input(base_inds(c,1):base_inds(c,2)),2);
            %disp(mean(input(base_inds(c,1):base_inds(c,2)),2));
            
            %pca
            spike_data = input(base_inds(c,1):base_inds(c,2));
            spike_data = spike_data - repmat(mean(spike_data,2),1,size(spike_data,2));
            [W, EvalueMatrix] = eig(cov(spike_data'));
            Evalues = diag(EvalueMatrix);
            Evalues = Evalues(end:-1:1);
            W = W(:,end:-1:1); W=W'; 
            pc = W * spike_data;
            %plot(pc(1,:),pc(2,:)) 
            %figure(4);
            %disp(pc(1));
            descending_indexes_of_local_maxs_bases(c,5) = pc(1);
            descending_indexes_of_local_maxs_bases(c,6) = pc(2);
            
            
            %remove submat of baseinds
            input(base_inds(c,1):base_inds(c,2)) = 0;
            
            c=c+1;
        end
     
    end
        
end


% function input = input_file_read(x)
% end

function  value = extract_features_cluster(rdf, der_input)
    %plot max vs var.
        cluster_coeff = 1;
        %figure('Color',[1 1 1]);
        figure(5);
        scatter(der_input(rdf(:,1)), rdf(:,4) );
        %scatter(der_input(rdf(:,1)), rdf(:,4) ,'k*');
         kmeans(rdf,cluster_coeff);
        hold all
        str= sprintf('max vs. var'); xlabel('max'); ylabel('var');title(str);%xlim([time(1) time(size(time,1))]);
        
%         %plot max vs var. unfiltered
%         figure(3);
%         scatter(test_input(rdf(:,1)), rdf(:,4) );
%         %scatter(der_input(rdf(:,1)), rdf(:,4) ,'k*');
%         kmeans(rdf,cluster_coeff);
%         hold all
%         str= sprintf('max vs. var (unfiltered)'); xlabel('max');title(str);%xlim([time(1) time(size(time,1))]);
 
    %plot pca1 vs pca 2
        cluster_coeff = 1;figure(6); %figure('Color',[1 1 1]);
        scatter(rdf(:,5), rdf(:,6)); kmeans(rdf,cluster_coeff);
        hold all
        str= sprintf('pca components'); xlabel('pca1'); ylabel('pca2');title(str);
    
     %plot max vs width
        cluster_coeff = 1; figure(7); %figure('Color',[1 1 1]);  
        scatter(der_input(rdf(:,1)), rdf(:,3) - rdf(:,2)); kmeans(rdf,cluster_coeff);
        hold all
        str= sprintf('max vs. width'); xlabel('spike max'); ylabel('spike width');title(str);%xlim([time(1) time(size(time,1))]); 
     
    %plot mean vs max
        cluster_coeff = 1;  figure(8); %figure('Color',[1 1 1]); 
        scatter(rdf(:,7), der_input(rdf(:,1) )); kmeans(rdf,cluster_coeff);
        hold all
        str= sprintf('mean vs max'); xlabel('mean'); ylabel('max');title(str);%xlim([time(1) time(size(time,1))]); 
        
        value = 1; %return value
end


function test_input = input_read(x)
    M= csvread(x); %read in csv file
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
end




% function filtered_input = filter_signal(input)
%     %derivative analysis
%     der_col = zeros(1,numel(input));
%     for idx= 2:numel(input)
%         der_col(idx) = input(idx)-input(idx-1); 
%     end
%     filtered_input = der_col;
% %     clc;
% %     plot(first_col, der_col);
% %     axis([0,21,-1,1]);
% end




%     base_inds(1,1) = min( abs(loc_max_ind - bases ));
%     dist_to_nearest_zero = min(abs(loc_max_ind - bases));
%     if( input (loc_max_ind - dist_to_nearest_zero) == 0)
%         base_inds(1,1) = loc_max_ind- dist_to_nearest_zero;
%     end
%     if (input (loc_max_ind + dist_to_nearest_zero) == 0)
%          base_inds(1,1) = loc_max_ind + dist_to_nearest_zero;
%     end
%     if( (input (loc_max_ind + dist_to_nearest_zero) == 0) && (input (max - dist_to_nearest_zero) == 0) )
%         disp('equidistant bases');
%         base_inds(1,1) = loc_max_ind + dist_to_nearest_zero;
%         base_inds(1,2) = loc_max_ind - dist_to_nearest_zero;
%     end
%        