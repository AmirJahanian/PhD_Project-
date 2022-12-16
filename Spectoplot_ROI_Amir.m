% ************************************************************************
% The code was developed at the Neurosceince Lab, Constructor University (Jacobs University)

% Contributors to MADE pipeline:
% Amir Jahanian Najafabadi (a.jahaniannajafabadi@jacobs-university.de)

% Please cite the following references for in any manuscripts produced utilizing MADE code:

% EEGLAB: A Delorme & S Makeig (2004) EEGLAB: an open source toolbox for
% analysis of single-trial EEG dynamics. Journal of Neuroscience Methods, 134, 9?21.

% This code is released under the GNU General Public License version 3.

% ************************************************************************

clear % clear matlab workspace
clc % clear matlab command window

% Enter the path of the file that has the data to be analyzed
rawdata_location = 'path/to/data/file.csv';

% Enter the path of the file where you want to save the processed data
output_location = 'path/to/output/file.csv';
output_dir = 'path/to/output/diractory/';
plot_location = 'path/to/plot/diractory';

channel_locations = 'path/to/channel/location/file';

ROI = [["F3" ,"F4", "Fz"], ["C3" , "C4", "Cz"], ["P3", "Pz", "P4"],["O1", "Oz", "O2"]];
nROI = 4;
ch_in_ROI = 3
freqfac = 10;
overlap = 0.5;
freqRange = [1 40];
cols_names = {"Frequency", "F", "C", "P", "O"};

sampling_rate = 256;
epoch_rate = 110;
ch_number = 12;

%% Read files to analyses
datafile_names=dir([rawdata_location filesep '*.mat']);
datafile_names=datafile_names(~ismember({datafile_names.name},{'.', '..', '.DS_Store'}));
datafile_names={datafile_names.name};
[filepath,name,ext] = fileparts(char(datafile_names{1}));

%%
AVG_EEG = zeros(ch_number, sampling_rate, epoch_rate)

%% Loop over all data files
for subject=1:length(datafile_names)
    
    EEG=[];
    ALLEEG = [];
    
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', subject, datafile_names{subject});
    
    
    EEG = pop_importdata('dataformat','matlab','nbchan',0,'data',[rawdata_location filesep datafile_names{subject}],'srate',sampling_rate...
       ,'pnts',0,'xmin',0,'chanlocs',channel_locations);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 

    %% Import channel locations
    EEG = pop_chanedit(EEG, 'lookup',channel_locations);
    EEG = eeg_checkset( EEG );
    
    for ch_index = 1:ch_number
        ch_name = AVG_EEG.data.chanlocs(ch_index).labels;

        %%
        % make sure we are averaging th correct channels together
        matching_index = -1;
        for ch_index2=1:EEG.data.nbchan
            if ch_name == EEG.data.chanlocs(ch_index2).labels
                matching_index = ch_index2;
                break;
            end
        end
        if matching_index == -1
            continue;
        end
        
        for sample_index = 1: sampling_rate
            for epoch_index=1:epoch_rate
                AVG_EEG(ch_index,sample_index,epoch_index) = AVG_EEG(ch_index,sample_index,epoch_index) + ...
                                                        (EEG.data.data(matching_index,sample_index,epoch_index)/length(datafile_names));
           end
        end
    end
end % end of subject loop

%%
PSD_data_length = (sampling_rate/2)*freqfac+1
final_avg = zeros(4,PSD_data_length);

for i=1:ch_number
    [spectra,freqs] = spectopo(AVG_EEG(i,:,:), 0, sampling_rate,'freqfac', freqfac, 'freqrange', [freqRange(1), freqRange(2)], 'overlap', ceil(sampling_rate/overlap));
      
    ch_name = AVG_EEG.data.chanlocs(i).labels;

    %%
    % make sure we are averaging th correct channels together
    insert_index = 0;
    for j=1:ch_number
        if ROI(j) == ch_name
            insert_index = ceil(j/ch_in_ROI);
            break;
        end
    end

    for j = 1: PSD_data_length
        final_avg(insert_index,j) = final_avg(insert_index,j)+(spectra(j));
    end
      
end

for i=1:nROI
    for j = 1:PSD_data_length
        final_avg(i,j) = final_avg(i,j)/ch_in_ROI;
    end
end

%%
for i =1:PSD_data_length
    if freqs(i) >= freqRange(1)
        break;
    end
    Frequency(i) = freqs(i);
    F(i) = final_avg(1,i);
    C(i) = final_avg(2,i);
    P(i) = final_avg(3,i);
    O(i) = final_avg(4,i);
end

%%
T = table(transpose(Frequency), transpose(F),transpose(C),transpose(P),transpose(O),'VariableNames', cols_names);
writetable(T,output_location);



