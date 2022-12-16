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
% Age   Session    Frequency   F   C   P   O

clear % clear matlab workspace
clc % clear matlab command window

% Enter the path of the file that has the data to be analyzed
filename = 'path/to/data/file.csv';

% Enter the path of the file where you want to save the processed data
output_location = 'path/to/output/file.csv';
output_dir = 'path/to/output/diractory/';

%%
ROI = ["O","P","F","C"];
Age_Group = ["YA","OA"];
Session = ["Pre-test","Mid-test"];

index = 1;
for age=1:length(Age_Group)
    for sess=1:length(Session)
        data(index).Age_Group = Age_Group(age);
        data(index).Session = Session(sess);
        data(index).Odata = [];
        data(index).Pdata = [];
        data(index).Fdata = [];
        data(index).Cdata = [];
        data(index).freqs = [];
        index = index + 1;
    end
end

%%
T = readtable(filename,'NumHeaderdata',1);

%%
for line=1:length(data)
    data(line).freqs = table2array(T(table2array(T(:,1)) == data(line).Age_Group & table2array(T(:,2)) == data(line).Session,3));
    data(line).Odata = table2array(T(table2array(T(:,1)) == data(line).Age_Group & table2array(T(:,2)) == data(line).Session,7));
    data(line).Pdata = table2array(T(table2array(T(:,1)) == data(line).Age_Group & table2array(T(:,2)) == data(line).Session,6));
    data(line).Cdata = table2array(T(table2array(T(:,1)) == data(line).Age_Group & table2array(T(:,2)) == data(line).Session,5));
    data(line).Fdata = table2array(T(table2array(T(:,1)) == data(line).Age_Group  & table2array(T(:,2)) == data(line).Session,4));
end

%%

fig = figure;
hold on 
title('Spectrogram ROI: O')
xlabel('Frequency (Hz)') 
ylabel('Log power Spectral Density 10*Log_{10}(uV^2/Hz)') 
index = 1;
for line=1:length(data)
    color = 'blue';
    if data(line).Age_Group == "YA"
        color = 'red';
    end
    if data(line).Session == "Pre-test"
        plot(data(line).freqs,data(line).Odata, '--', 'color',color);
    else
        plot(data(line).freqs,data(line).Odata,  'color',color);
    end
    Okeys(index) = strjoin([data(line).Age_Group "_" data(line).Session]);
    index = index + 1;
end
xlim([0 40]);
ylim([-25 25]);
legend(Okeys);
hold off
filename = strjoin([output_dir 'PSD_plot_roi_O.png']);
saveas(fig,filename);

fig = figure;
hold on 
title('Spectrogram ROI: P')
xlabel('Frequency (Hz)') 
ylabel('Log power Spectral Density 10*Log_{10}(uV^2/Hz)') 
index = 1;
for line=1:length(data)
    color = 'blue';
    if data(line).Age_Group == "YA"
        color = 'red';
    end
    if data(line).Session == "Pre-test"
        plot(data(line).freqs,data(line).Pdata, '--', 'color',color);
    else
        plot(data(line).freqs,data(line).Pdata,  'color',color);
    end
    Pkeys(index) = strjoin([data(line).Age_Group "_" data(line).Session]);
    index = index + 1;
end
xlim([0 40]);
ylim([-25 25]);
legend(Pkeys);
hold off
filename = strjoin([output_dir 'PSD_plot_roi_P.png']);
saveas(fig,filename);

fig = figure;
hold on 
title('Spectrogram ROI: C')
xlabel('Frequency (Hz)') 
ylabel('Log power Spectral Density 10*Log_{10}(uV^2/Hz)') 
index = 1;
for line=1:length(data)
    color = 'blue';
    if data(line).Age_Group == "YA"
        color = 'red';
    end
    if data(line).Session == "Pre-test"
        plot(data(line).freqs,data(line).Cdata, '--', 'color',color);
    else
        plot(data(line).freqs,data(line).Cdata,  'color',color);
    end
    Ckeys(index) = strjoin([data(line).Age_Group "_" data(line).Session]);
    index = index + 1;
end
xlim([0 40]);
ylim([-25 25]);
legend(Ckeys);
hold off
filename = strjoin([output_dir 'PSD_plot_roi_C.png']);
saveas(fig,filename);

fig = figure;
hold on 
title('Spectrogram ROI: F')
xlabel('Frequency (Hz)') 
ylabel('Log power Spectral Density 10*Log_{10}(uV^2/Hz)') 
index = 1;
for line=1:length(data)
    color = 'blue';
    if data(line).Age_Group == "YA"
        color = 'red';
    end
    if data(line).Session == "Pre-test"
        plot(data(line).freqs,data(line).Fdata, '--', 'color',color);
    else
        plot(data(line).freqs,data(line).Fdata,  'color',color);
    end
    Fkeys(index) = strjoin([data(line).Age_Group "_" data(line).Session]);
    index = index + 1;
end
xlim([0 40]);
ylim([-25 25]);
legend(Fkeys);
hold off
filename = strjoin([output_dir 'PSD_plot_roi_F.png']);
saveas(fig,filename);

