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


output_location = 'C:\PALMS Elderly REEG\Spectogram_Processed_datasets\OA\first_AR\V\processed_data';

filename_task_V = ['L:\PALMS EEG Processing' filesep 'OA_EEG_AR_V_Finalversion.xls'];
filename_task_VT = ['L:\PALMS EEG Processing' filesep 'OA_EEG_AR_VT_Finalversion.xls'];

filename_REEG = ['D:\Amir_YA_OA_MATLAB_Code' filesep 'OA_REEG_Pre-test_Finalversion35.xls'];

ROI = ["O","P","F","C"];

Age_Group = ["YA","OA"];

Feedback_Type = ["V","VT"];
Session = ["Pre-test"];

index = 1;
for ft=1:length(Feedback_Type)
    lines(index).Feedback_Type = Feedback_Type(ft);
    lines(index).Odata = [];
    lines(index).Pdata = [];
    lines(index).Fdata = [];
    lines(index).Cdata = [];
    lines(index).freqs = [];
    index = index + 1;
        
end

for ag=1:length(Age_Group)
    lines_REEG(index).Age_Group = Age_Group(ag);
    lines_REEG(index).Odata = [];
    lines_REEG(index).Pdata = [];
    lines_REEG(index).Fdata = [];
    lines_REEG(index).Cdata = [];
    lines_REEG(index).freqs = [];
    index = index + 1;
end


%%
T_V = readtable(filename_task_V,'NumHeaderLines',1);

T_VT = readtable(filename_task_VT,'NumHeaderLines',1);
T_REEG = readtable(filename_REEG,'NumHeaderLines',1);


%%

% Age	Session 	Feedback	Frequency	F	C	P	O


for line=1:length(lines)
    if lines(line).Feedback_Type == "VT"
        lines(line).freqs = table2array(T_VT(:,1));
        lines(line).Odata = table2array(T_VT(:,5));
        lines(line).Pdata = table2array(T_VT(:,4));
        lines(line).Cdata = table2array(T_VT(:,3));
        lines(line).Fdata = table2array(T_VT(:,2));
    else
        lines(line).freqs = table2array(T_V(:,1));
        lines(line).Odata = table2array(T_V(:,5));
        lines(line).Pdata = table2array(T_V(:,4));
        lines(line).Cdata = table2array(T_V(:,3));
        lines(line).Fdata = table2array(T_V(:,2));
    end
    
    
end


for line=1:length(lines_REEG)
    lines_REEG(line).freqs = table2array(T_REEG(:,1));
    lines_REEG(line).Odata = table2array(T_REEG(:,5));
    lines_REEG(line).Pdata = table2array(T_REEG(:,4));
    lines_REEG(line).Cdata = table2array(T_REEG(:,3));
    lines_REEG(line).Fdata = table2array(T_REEG(:,2));    
    
end

%%
% Age	Session 	Feedback	Frequency	F	C	P	O

fig = figure;
hold on 
title('Spectrogram ROI: O')
xlabel('Frequency (Hz)') 
ylabel('Log power Spectral Density 10*Log_{10}(uV^2/Hz)') 
index = 1;
k=1;
for line=1:length(lines)
        if index == 1
            plot(lines_REEG(k).freqs,lines_REEG(k).Odata,  'color','blue');
            Okeys(index) = "REEG";
            index = index + 1;
        end
            
        if lines(line).Feedback_Type == "V"
            plot(lines(line).freqs,minus(lines(line).Odata,lines_REEG(k).Odata), '--',  'color','red');
        else
            plot(lines(line).freqs,minus(lines(line).Odata,lines_REEG(k).Odata),  'color','red');
        end
            

        Okeys(index) = lines(line).Feedback_Type;
        index = index + 1;
end

xlim([0 30]);
ylim([-25 25]);
legend(Okeys);
hold off
filename = 'output_roi_O.png';
saveas(fig,filename);






fig = figure;
hold on 
title('Spectrogram ROI: P')
xlabel('Frequency (Hz)') 
ylabel('Log power Spectral Density 10*Log_{10}(uV^2/Hz)') 
index = 1;
for line=1:length(lines)
        if index == 1
            plot(lines_REEG(k).freqs,lines_REEG(k).Pdata,  'color','blue');
            Pkeys(index) = "REEG";
            index = index + 1;
        end
            
        if lines(line).Feedback_Type == "V"
            plot(lines(line).freqs,minus(lines(line).Pdata,lines_REEG(k).Pdata), '--',  'color','red');
        else
            plot(lines(line).freqs,minus(lines(line).Pdata,lines_REEG(k).Pdata),  'color','red');
        end
            

        Pkeys(index) = lines(line).Feedback_Type;
        index = index + 1;
end

xlim([0 30]);
ylim([-25 25]);
legend(Pkeys);
hold off
filename = 'output_roi_P.png';
saveas(fig,filename);








fig = figure;
hold on 
title('Spectrogram ROI: C')
xlabel('Frequency (Hz)') 
ylabel('Log power Spectral Density 10*Log_{10}(uV^2/Hz)') 
index = 1;
for line=1:length(lines)
        if index == 1
            plot(lines_REEG(k).freqs,lines_REEG(k).Cdata,  'color','blue');
            Ckeys(index) = "REEG";
            index = index + 1;
        end
            
        if lines(line).Feedback_Type == "V"
            plot(lines(line).freqs,minus(lines(line).Cdata,lines_REEG(k).Cdata), '--',  'color','red');
        else
            plot(lines(line).freqs,minus(lines(line).Cdata,lines_REEG(k).Cdata),  'color','red');
        end
            

        Ckeys(index) = lines(line).Feedback_Type;
        index = index + 1;
end

xlim([0 30]);
ylim([-25 25]);
legend(Ckeys);
hold off
filename = 'output_roi_C.png';
saveas(fig,filename);







fig = figure;
hold on 
title('Spectrogram ROI: F')
xlabel('Frequency (Hz)') 
ylabel('Log power Spectral Density 10*Log_{10}(uV^2/Hz)') 
index = 1;
for line=1:length(lines)
        if index == 1
            plot(lines_REEG(k).freqs,lines_REEG(k).Fdata,  'color','blue');
            Fkeys(index) = "REEG";
            index = index + 1;
        end
            
        if lines(line).Feedback_Type == "V"
            plot(lines(line).freqs,minus(lines(line).Fdata,lines_REEG(k).Fdata), '--',  'color','red');
        else
            plot(lines(line).freqs,minus(lines(line).Fdata,lines_REEG(k).Fdata),  'color','red');
        end
            

        Fkeys(index) = lines(line).Feedback_Type;
        index = index + 1;
end

xlim([0 30]);
ylim([-25 25]);
legend(Fkeys);
hold off
filename = 'OA_roi_F_REEG_Feedback.png';
saveas(fig,filename);





