% ************************************************************************
% The code was developed at the Neurosceince Lab, Constructor University (Jacobs University)

% Contributors to MADE pipeline:
% Amir Jahanian Najafabadi (a.jahaniannajafabadi@jacobs-university.de)

% Please cite the following references for in any manuscripts produced utilizing MADE code:

% EEGLAB: A Delorme & S Makeig (2004) EEGLAB: an open source toolbox for
% analysis of single-trial EEG dynamics. Journal of Neuroscience Methods, 134, 9?21.

% This code is released under the GNU General Public License version 3.

% ************************************************************************

% The following code reads processed PSD values that were averaged acorrding to ROI 
% and plot the results together according to ROI, age group, and session
% The assumption on the input file as as follows:
% The file is either: csv,xls,xlsx
% The file has the following columns structure
% Age   Session     Feedback    Frequency   F   C   P   O


clear % clear matlab workspace
clc % clear matlab command window

% Enter the path of the file that has the data to be analyzed
filename = 'path/to/data/file.csv';

% Enter the path of the file where you want to save the processed data
output_location = 'path/to/output/file.csv';
output_dir = 'path/to/output/diractory/';

ROI = ["O","P","F","C"];
Age_Group = ["YA","OA"];
Sessions = [1,2];
Feedback_Type = ["V","VT"];

index = 1;
for ag=1:length(Age_Group)
    for sess=1:length(Feedback_Type)
        for session=1:length(Sessions)
            data(index).Age_Group = Age_Group(ag);
            data(index).Feedback_Type = Feedback_Type(sess);
            data(index).Session = Sessions(session);
            data(index).Odata = [];
            data(index).Pdata = [];
            data(index).Fdata = [];
            data(index).Cdata = [];
            data(index).freqs = [];
            index = index + 1;
        end
    end
end

%%
T = readtable(filename,'NumHeaderdata',1);

%%
for line=1:length(data)
    data(line).freqs = table2array(T(  table2array(T(:,1)) == data(line).Age_Group & ...
                                        table2array(T(:,2)) == data(line).Session & ...
                                        table2array(T(:,3)) == data(line).Feedback_Type,4));

    data(line).Odata = table2array(T(  table2array(T(:,1)) == data(line).Age_Group & ...
                                        table2array(T(:,2)) == data(line).Session & ...
                                        table2array(T(:,3)) == data(line).Feedback_Type,8));

    data(line).Pdata = table2array(T(  table2array(T(:,1)) == data(line).Age_Group & ...
                                        table2array(T(:,2)) == data(line).Session & ...
                                        table2array(T(:,3)) == data(line).Feedback_Type,7));

    data(line).Cdata = table2array(T(  table2array(T(:,1)) == data(line).Age_Group & ...
                                        table2array(T(:,2)) == data(line).Session & ...
                                        table2array(T(:,3)) == data(line).Feedback_Type,6));

    data(line).Fdata = table2array(T(  table2array(T(:,1)) == data(line).Age_Group & ...
                                        table2array(T(:,2)) == data(line).Session & ...
                                        table2array(T(:,3)) == data(line).Feedback_Type,5));
end

%%
for session=1:length(Sessions)
    fig = figure;
    hold on 
    title(['Spectrogram ROI: O, Session: ' string(session)])
    xlabel('Frequency (Hz)') 
    ylabel('Log power Spectral Density 10*Log_{10}(uV^2/Hz)') 
    index = 1;
    for line=1:length(data)
        if data(line).Session == session
            color = 'blue';
            if data(line).Age_Group == "YA"
                color = 'red';
            end
            if data(line).Feedback_Type == "V"
                plot(data(line).freqs,data(line).Odata, '--', 'color',color);
            else
                plot(data(line).freqs,data(line).Odata,  'color',color);
            end
            Okeys(index) = strjoin([data(line).Age_Group "_" data(line).Feedback_Type]);
            index = index + 1;
        end
    end
    xlim([0 40]);
    ylim([-25 25]);
    legend(Okeys);
    hold off
    filename = strjoin([output_dir 'PSD_plot_roi_O_session_' string(session) '.png']);
    saveas(fig,filename);

    fig = figure;
    hold on 
    title(['Spectrogram ROI: P, Session: ' string(session)])
    xlabel('Frequency (Hz)') 
    ylabel('Log power Spectral Density 10*Log_{10}(uV^2/Hz)') 
    index = 1;
    for line=1:length(data)
        if data(line).Session == session
            color = 'blue';
            if data(line).Age_Group == "YA"
                color = 'red';
            end
            if data(line).Feedback_Type == "V"
                plot(data(line).freqs,data(line).Pdata, '--', 'color',color);
            else
                plot(data(line).freqs,data(line).Pdata,  'color',color);
            end
            Pkeys(index) = strjoin([data(line).Age_Group "_" data(line).Feedback_Type]);
            index = index + 1;
        end
    end
    xlim([0 40]);
    ylim([-25 25]);
    legend(Pkeys);
    hold off
    filename = strjoin([output_dir 'PSD_plot_roi_P_session_' string(session) '.png']);
    saveas(fig,filename);

    fig = figure;
    hold on 
    title(['Spectrogram ROI: C, Session: ' string(session)])
    xlabel('Frequency (Hz)') 
    ylabel('Log power Spectral Density 10*Log_{10}(uV^2/Hz)') 
    index = 1;
    for line=1:length(data)
        if data(line).Session == session
            color = 'blue';
            if data(line).Age_Group == "YA"
                color = 'red';
            end
            if data(line).Feedback_Type == "V"
                plot(data(line).freqs,data(line).Cdata, '--', 'color',color);
            else
                plot(data(line).freqs,data(line).Cdata,  'color',color);
            end
            Ckeys(index) = strjoin([data(line).Age_Group "_" data(line).Feedback_Type]);
            index = index + 1;
        end
    end
    xlim([0 40]);
    ylim([-25 25]);
    legend(Ckeys);
    hold off
    filename = strjoin([output_dir 'PSD_plot_roi_C_session_' string(session) '.png']);
    saveas(fig,filename);

    fig = figure;
    hold on 
    title(['Spectrogram ROI: F, Session: ' string(session)])
    xlabel('Frequency (Hz)') 
    ylabel('Log power Spectral Density 10*Log_{10}(uV^2/Hz)') 
    index = 1;
    for line=1:length(data)
        if data(line).Session == session
            color = 'blue';
            if data(line).Age_Group == "YA"
                color = 'red';
            end
            if data(line).Feedback_Type == "V"
                plot(data(line).freqs,data(line).Fdata, '--', 'color',color);
            else
                plot(data(line).freqs,data(line).Fdata,  'color',color);
            end
            Fkeys(index) = strjoin([data(line).Age_Group "_" data(line).Feedback_Type]);
            index = index + 1;
        end
    end
    xlim([0 40]);
    ylim([-25 25]);
    legend(Fkeys);
    hold off
    filename = strjoin([output_dir 'PSD_plot_roi_O_session_' string(session) '.png']);
    saveas(fig,filename);
end
