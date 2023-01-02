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

% data path
rawdata_location_V = 'path/to/data/with/feedback/V';
rawdata_location_VT = 'path/to/data/with/feedback/VT';
rawdata_location_Posttest = 'path/to/data/REEG/pretest';
rawdata_location_Pretest = 'path/to/data/REEG/posttest';
channel_locations = 'path/to/channel/location';

% data parameters
sampling_rate = 256;
epoch_rate = 120;
ch_number = 12;

channel_name = {'F3', 'C3', 'P3', 'Pz', 'O1', 'Oz', 'O2', 'P4', 'C4', 'F4', 'Fz', 'Cz'};

% this is a list of participant ids that had task related V in training in 
% their first session (corresponding to pre test REEG)
% this parameter is used to generate the correct averaging between task ralted and REEG
Id_map_OA_V_preTest = ["PEO8","PEO10","PEO12","PEO14","PEO16","PEO18","PEO20","PEO22","PEO24","PEO27","PEO28","PEO32","PEO34","PEO35","PEO38","PEO39","PEO41"];


% spectotopo parameters
freqfac = 10;
freqrange = [1,40];
overlap = ceil(sampling_rate/2);

% ROI parameters (do not change these parameters, there are hard coded 
% parts of the code dependent on these values)
nROI = 4;
chennel_per_ROI = 3; 
ROI = [["F3" ,"F4", "Fz"], ["C3" , "C4", "Cz"], ["P3", "Pz", "P4"],["O1", "Oz", "O2"]];

%%

datafile_names_V=dir([rawdata_location_V filesep '*.mat']);
datafile_names_V=datafile_names_V(~ismember({datafile_names_V.name},{'.', '..', '.DS_Store'}));
datafile_names_V={datafile_names_V.name};
[filepath,name,ext] = fileparts(char(datafile_names_V{1}));



datafile_names_VT=dir([rawdata_location_VT filesep '*.mat']);
datafile_names_VT=datafile_names_VT(~ismember({datafile_names_VT.name},{'.', '..', '.DS_Store'}));
datafile_names_VT={datafile_names_VT.name};
[filepath,name,ext] = fileparts(char(datafile_names_VT{1}));



datafile_names_Pre=dir([rawdata_location_Pretest filesep '*.mat']);
datafile_names_Pre=datafile_names_Pre(~ismember({datafile_names_Pre.name},{'.', '..', '.DS_Store'}));
datafile_names_Pre={datafile_names_Pre.name};
[filepath,name,ext] = fileparts(char(datafile_names_Pre{1}));


datafile_names_Post=dir([rawdata_location_Posttest filesep '*.mat']);
datafile_names_Post=datafile_names_Post(~ismember({datafile_names_Post.name},{'.', '..', '.DS_Store'}));
datafile_names_Post={datafile_names_Post.name};
[filepath,name,ext] = fileparts(char(datafile_names_Post{1}));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop over data with feedback v, compute the spectogram and place all the results in struct V_EEG
index = 1
for subject=1:length(datafile_names_V)

    EEG=[];
    ALLEEG = [];

    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', subject, datafile_names_V{subject});

    EEG = pop_importdata('dataformat','matlab','nbchan',0,'data',[rawdata_location_V filesep datafile_names_V{subject}],'srate',sampling_rate...
       ,'pnts',0,'xmin',0,'chanlocs',channel_locations);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off');

    EEG=pop_chanedit(EEG, 'lookup',channel_locations);
    EEG = eeg_checkset( EEG );
    [spectra freqs] = spectopo(EEG.data.data, 0, sampling_rate,'freqfac', freqfac, 'freqrange', freqrange, 'overlap', overlap);
    V_EEG(index).data = spectra
    V_EEG(index).freqs = freqs;
    V_EEG(index).n_ch = EEG.data.nbchan;
    for ch=1:length(EEG.data.chanlocs)
        V_EEG(index).ch_name(ch) = string(EEG.data.chanlocs(ch).labels);
    end
    id = split(datafile_names_V(subject),"_");
    V_EEG(index).id = string(id(1));
    index = index + 1;
end % end of subject loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop over data with feedback vt, compute the spectogram and place all the results in struct VT_EEG
index = 1
for subject=1:length(datafile_names_VT)

    EEG=[];
    ALLEEG = [];

    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', subject, datafile_names_VT{subject});

    EEG = pop_importdata('dataformat','matlab','nbchan',0,'data',[rawdata_location_VT filesep datafile_names_VT{subject}],'srate',sampling_rate...
       ,'pnts',0,'xmin',0,'chanlocs',channel_locations);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off');

    EEG=pop_chanedit(EEG, 'lookup',channel_locations);
    EEG = eeg_checkset( EEG );
    [spectra freqs] = spectopo(EEG.data.data, 0, sampling_rate,'freqfac', freqfac, 'freqrange', freqrange, 'overlap', overlap);
    VT_EEG(index).data = spectra
    VT_EEG(index).freqs = freqs;
    VT_EEG(index).n_ch = EEG.data.nbchan;
    for ch=1:length(EEG.data.chanlocs)
        VT_EEG(index).ch_name(ch) = string(EEG.data.chanlocs(ch).labels);
    end
    id = split(datafile_names_VT(subject),"_");
    VT_EEG(index).id = string(id(1));
    index = index + 1;
end % end of subject loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop over data with REEG pretest, compute the spectogram and place all the results in struct Pre_EEG
index = 1
for subject=1:length(datafile_names_Pre)

    EEG=[];
    ALLEEG = [];

    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', subject, datafile_names_Pre{subject});

    EEG = pop_importdata('dataformat','matlab','nbchan',0,'data',[rawdata_location_Pretest filesep datafile_names_Pre{subject}],'srate',sampling_rate...
       ,'pnts',0,'xmin',0,'chanlocs',channel_locations);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off');

    EEG=pop_chanedit(EEG, 'lookup',channel_locations);
    EEG = eeg_checkset( EEG );
    [spectra freqs] = spectopo(EEG.data.data, 0, sampling_rate,'freqfac', freqfac, 'freqrange', freqrange, 'overlap', overlap);
    Pre_EEG(index).data = spectra
    Pre_EEG(index).freqs = freqs;
    Pre_EEG(index).n_ch = EEG.data.nbchan;
    for ch=1:length(EEG.data.chanlocs)
        Pre_EEG(index).ch_name(ch) = string(EEG.data.chanlocs(ch).labels);
    end
    id = split(datafile_names_Pre(subject),"_");
    Pre_EEG(index).id = string(id(1));
    index = index + 1;
end % end of subject loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop over data with REEG posttest, compute the spectogram and place all the results in struct Post_EEG
index = 1
for subject=1:length(datafile_names_Post)

    EEG=[];
    ALLEEG = [];

    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', subject, datafile_names_Post{subject});

    EEG = pop_importdata('dataformat','matlab','nbchan',0,'data',[rawdata_location_Posttest filesep datafile_names_Post{subject}],'srate',sampling_rate...
       ,'pnts',0,'xmin',0,'chanlocs',channel_locations);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off');

    EEG=pop_chanedit(EEG, 'lookup',channel_locations);
    EEG = eeg_checkset( EEG );
    [spectra freqs] = spectopo(EEG.data.data, 0, sampling_rate,'freqfac', freqfac, 'freqrange', freqrange, 'overlap', overlap);
    Post_EEG(index).data = spectra
    Post_EEG(index).freqs = freqs;
    Post_EEG(index).n_ch = EEG.data.nbchan;
    for ch=1:length(EEG.data.chanlocs)
        Post_EEG(index).ch_name(ch) = string(EEG.data.chanlocs(ch).labels);
    end
    id = split(datafile_names_Post(subject),"_");
    Post_EEG(index).id = string(id(1));
    index = index + 1;
end % end of subject loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% number of samples outputed by spectotopo function = (sampling_rate/2)*freqfac+1

% note: for REEG averaging, both post and pre test are averaged together for the same frequency
REEG = zeros(ch_number, ((sampling_rate/2)*freqfac+1) );
TR_V = zeros(ch_number, ((sampling_rate/2)*freqfac+1) );
TR_VT = zeros(ch_number, ((sampling_rate/2)*freqfac+1) );

%% Loop over V_EEG and use the map provided in parameters to compute the correct averaging
for v=1:length(V_EEG)
    % look for the currect id in the ma provided
    found_id = false;
    %% REEG
    for map=1:length(Id_map_OA_V_preTest)
       if string(V_EEG(v).id(1) )== Id_map_OA_V_preTest(map)
            found_id = true;
            break
        end
    end
    reeg = struct
    add = false;
    % if found, then subtract pretest from the current sample
    if found_id
        for REEG_subject=1:length(Pre_EEG)
            if V_EEG(v).id(1) == Pre_EEG(REEG_subject).id(1)
                reeg = Pre_EEG(REEG_subject);
                add = true;
                break
            end
        end
    else % else subtract posttest from current sample
        for REEG_subject=1:length(Post_EEG)
            if V_EEG(v).id(1) == Post_EEG(REEG_subject).id(1)
                reeg = Post_EEG(REEG_subject);
                add = true;
                break
            end
        end
    end
    % if no match was found between task related and REEG, then skip
    if ~add 
        continue
    end

    % loop over channels and make sure to subtract the correct task 
    % related channel from REEG channel
    for ch_index = 1:length(channel_name)
        ch_name = string(channel_name(ch_index));
        %%
        matching_index = -1;
        for ch_index2=1:V_EEG(v).n_ch
            if ch_name == V_EEG(v).ch_name(ch_index2)
                matching_index = ch_index2;
                break;
            end
        end

        if matching_index == -1
            continue;
        end

        %%
        REEGmatching_index = -1;
        for ch_index2=1:reeg.n_ch
            if ch_name == reeg.ch_name(ch_index2)
                REEGmatching_index = ch_index2;
                break;
            end
        end

        %%

        for j = 1: ((sampling_rate/2)*freqfac+1)              
            datapoint = V_EEG(v).data(matching_index,j)-reeg.data(REEGmatching_index,j);
            TR_V(ch_index,j) = TR_V(ch_index,j) + (datapoint/ length(V_EEG));
        end
    end

end


%% Loop over VT_EEG and use the map provided in parameters to compute the correct averaging
for v=1:length(VT_EEG)
    found_id = false;
    % look for the currect id in the map provided
    for map=1:length(Id_map_OA_V_preTest)
       if VT_EEG(v).id(1) == Id_map_OA_V_preTest(map)
            found_id = true;
            break
        end
    end
    add = false;
    reeg = struct
    % if found, then subtract pretest from the current sample
    if ~found_id
        for REEG_subject=1:length(Pre_EEG)
            if VT_EEG(v).id(1)  == Pre_EEG(REEG_subject).id(1) 
                reeg = Pre_EEG(REEG_subject);
                add = true;
                break
            end
        end
    else % else subtract posttest from current sample
        for REEG_subject=1:length(Post_EEG)
            if VT_EEG(v).id(1)  == Post_EEG(REEG_subject).id(1) 
                reeg = Post_EEG(REEG_subject);
                add = true;
                break
            end
        end
    end
    % if no match was found between task related and REEG, then skip
    if ~add 
        continue
    end

    % loop over channels and make sure to subtract the correct task 
    % related channel from REEG channel
    for ch_index = 1:length(channel_name)
        ch_name = string(channel_name(ch_index));
        %%
        matching_index = -1;
        for ch_index2=1:VT_EEG(v).n_ch
            if ch_name == VT_EEG(v).ch_name(ch_index2)
                matching_index = ch_index2;
                break;
            end
        end

        if matching_index == -1
            continue;
        end

        %%
        REEGmatching_index = -1;
        for ch_index2=1:reeg.n_ch
            if ch_name == reeg.ch_name(ch_index2)
                REEGmatching_index = ch_index2;
                break;
            end
        end
        
        if REEGmatching_index == -1
            continue;
        end

        %%

        for j = 1: ((sampling_rate/2)*10+1)             
            datapoint = VT_EEG(v).data(matching_index,j)-reeg.data(REEGmatching_index,j);
            TR_VT(ch_index,j) = TR_VT(ch_index,j) + (datapoint/ length(VT_EEG));
        end
    end

end


%% Loop over Post_EEG and use the map provided in parameters to compute the correct averaging
for v=1:length(Post_EEG)
    for ch_index = 1:length(channel_name)
        ch_name = string(channel_name(ch_index));
        %%
        matching_index = -1;
        for ch_index2=1:Post_EEG(v).n_ch
            if ch_name == Post_EEG(v).ch_name(ch_index2)
                matching_index = ch_index2;
                break;
            end
        end

        if matching_index == -1
            continue;
        end

        %%

        for j = 1: ((sampling_rate/2)*freqfac+1)             
            datapoint = Post_EEG(v).data(matching_index,j);
            REEG(ch_index,j) = REEG(ch_index,j) + (datapoint/ (length(Post_EEG)+length(Pre_EEG)));
        end
    end

end

%% Loop over Pre_EEG and use the map provided in parameters to compute the correct averaging
for v=1:length(Pre_EEG)
    for ch_index = 1:length(channel_name)
        ch_name = string(channel_name(ch_index));
        %%
        matching_index = -1;
        for ch_index2=1:Pre_EEG(v).n_ch
            if ch_name == Pre_EEG(v).ch_name(ch_index2)
                matching_index = ch_index2;
                break;
            end
        end

        if matching_index == -1
            continue;
        end

        %%

        for j = 1: ((sampling_rate/2)*freqfac+1)             
            datapoint = Pre_EEG(v).data(matching_index,j);
            REEG(ch_index,j) = REEG(ch_index,j) + (datapoint/ (length(Pre_EEG)+length(Post_EEG)));
        end
    end

end


%%


final_avg_tr_V = zeros(nROI,(sampling_rate/2)*freqfac+1);
final_avg_tr_VT = zeros(nROI,(sampling_rate/2)*freqfac+1);
final_avg_reeg = zeros(nROI,(sampling_rate/2)*freqfac+1);

% compute final average
for i=1:ch_number
    ch_name = string(channel_name(i));
    insert_index = 0;
    for j=1:ch_number
        if ROI(j) == ch_name
            insert_index = ceil(j/chennel_per_ROI);
            break;
        end
    end
    for j = 1: ((sampling_rate/2)*freqfac+1)
        final_avg_tr_V(insert_index,j) = final_avg_tr_V(insert_index,j)+((TR_V(i,j))/chennel_per_ROI);
        final_avg_tr_VT(insert_index,j) = final_avg_tr_VT(insert_index,j)+((TR_VT(i,j))/chennel_per_ROI);
        final_avg_reeg(insert_index,j) = final_avg_reeg(insert_index,j)+((REEG(i,j))/chennel_per_ROI);
    end
end

%%
% only take data points up to 30 Hz (until the frequency of intrest)
for i =1:((sampling_rate/2)*freqfac+1)
    if freqs(i) >= 30
        break;
    end

    Frequency(i) = freqs(i);
    F_V(i) = final_avg_tr_V(1,i);
    C_V(i) = final_avg_tr_V(2,i);
    P_V(i) = final_avg_tr_V(3,i);
    O_V(i) = final_avg_tr_V(4,i);
end



for i =1:((sampling_rate/2)*freqfac+1)
    if freqs(i) >= 30
        break;
    end

    F_VT(i) = final_avg_tr_VT(1,i);
    C_VT(i) = final_avg_tr_VT(2,i);
    P_VT(i) = final_avg_tr_VT(3,i);
    O_VT(i) = final_avg_tr_VT(4,i);
end


for i =1:((sampling_rate/2)*freqfac+1)
    if freqs(i) >= 30
        break;
    end

    Frequency(i) = freqs(i);
    F(i) = final_avg_reeg(1,i);
    C(i) = final_avg_reeg(2,i);
    P(i) = final_avg_reeg(3,i);
    O(i) = final_avg_reeg(4,i);
end



%%
% save data into csv files
T_V = table(transpose(Frequency), transpose(F_V),transpose(C_V),transpose(P_V),transpose(O_V),'VariableNames', ["Frequency", "F", "C", "P", "O"]);

writetable(T_V,'V.xls');

T_VT = table(transpose(Frequency), transpose(F_VT),transpose(C_VT),transpose(P_VT),transpose(O_VT),'VariableNames', ["Frequency", "F", "C", "P", "O"]);

writetable(T_VT,'VT.xls');

T = table(transpose(Frequency), transpose(F),transpose(C),transpose(P),transpose(O),'VariableNames', ["Frequency", "F", "C", "P", "O"]);

writetable(T,'reeg.xls');

%%





% plot data and save plots

fig = figure;
hold on 
title('Spectrogram ROI: F')
xlabel('Frequency (Hz)') 
ylabel('Log power Spectral Density 10*Log_{10}(uV^2/Hz)') 
index = 1;
plot(Frequency,F,  'color','blue');
Fkeys(index) = "REEG";
index = index + 1;
plot(Frequency,F_V, '--',  'color','red');
Fkeys(index) = "V";
index = index + 1;

plot(Frequency,F_VT,  'color','red');


Fkeys(index) = "VT";
index = index + 1;

xlim([0 30]);
ylim([-15 15]);
legend(Fkeys);
hold off
filename = 'FINAL_ROI_F.png';
saveas(fig,filename);






fig = figure;
hold on 
title('Spectrogram ROI: P')
xlabel('Frequency (Hz)') 
ylabel('Log power Spectral Density 10*Log_{10}(uV^2/Hz)') 
index = 1;
plot(Frequency,P,  'color','blue');
Pkeys(index) = "REEG";
index = index + 1;
plot(Frequency,P_V, '--',  'color','red');
Pkeys(index) = "V";
index = index + 1;

plot(Frequency,P_VT,  'color','red');


Pkeys(index) = "VT";
index = index + 1;

xlim([0 30]);
ylim([-15 15]);
legend(Pkeys);
hold off
filename = 'FINAL_ROI_P.png';
saveas(fig,filename);






fig = figure;
hold on 
title('Spectrogram ROI: O')
xlabel('Frequency (Hz)') 
ylabel('Log power Spectral Density 10*Log_{10}(uV^2/Hz)') 
index = 1;
plot(Frequency,O,  'color','blue');
Okeys(index) = "REEG";
index = index + 1;
plot(Frequency,O_V, '--',  'color','red');
Okeys(index) = "V";
index = index + 1;

plot(Frequency,O_VT,  'color','red');


Okeys(index) = "VT";
index = index + 1;

xlim([0 30]);
ylim([-15 15]);
legend(Okeys);
hold off
filename = 'FINAL_ROI_O.png';
saveas(fig,filename);






fig = figure;
hold on 
title('Spectrogram ROI: C')
xlabel('Frequency (Hz)') 
ylabel('Log power Spectral Density 10*Log_{10}(uV^2/Hz)') 
index = 1;
plot(Frequency,C,  'color','blue');
Ckeys(index) = "REEG";
index = index + 1;
plot(Frequency,C_V, '--',  'color','red');
Ckeys(index) = "V";
index = index + 1;

plot(Frequency,C_VT,  'color','red');


Ckeys(index) = "VT";
index = index + 1;

xlim([0 30]);
ylim([-15 15]);
legend(Fkeys);
hold off
filename = 'FINAL_ROI_C.png';
saveas(fig,filename);



















