%get LFP phase for all cells in the dataset
%ignore Cluster0 for all
%clear; close all; clc;
function getLFPphase(fID)

drMAT = '/scratch/dd1348/fraser/MAT/';
drTS = '/scratch/dd1348/fraser/TETRODE_TS/';
drART = '/scratch/dd1348/fraser/ART_border025dif/';
drOUT = '/scratch/dd1348/fraser/PHASES/';

%list of files to process
%fid = fopen('filesToProcess_10_04_2015.txt'); %28
%fid = fopen('filesToProcess_10_10_2015.txt'); %18
%fid = fopen('filesOpenField_1_31_2016.txt'); %21
%files = textscan(fid, '%s');
%fclose(fid);
%files = files{1};

files = dir([drMAT '*.mat']);
fn = files(fID).name;
disp(fn);

%skip if file exists
if exist([drOUT fn],'file') ~= 0; return; end;

%create wavelets for all bands
bands = logspace(log10(2)/log10(10),log10(250)/log10(10),100);
%bands = [20,30];
Nbands = length(bands);
eegFS = 2000;
wFactor = 8;
psiBands = getWavelets( bands, wFactor, eegFS );

%load files
load([drMAT fn]);
load([drTS fn]);
load([drART fn]);

%make timestamps for EEG
tsEEG = nan(1,size(eegData,2));
eegTimeStamp = eegTimeStamp/10; %ms
eegTimeStamp = double(eegTimeStamp);
Nblocks = size(eegData,2) / length(eegTimeStamp);
dT = 1000/eegFS; %sample time
ts = 1:Nblocks; %times within each block
ts = ts-1; %starting from zero
ts = ts*dT; %convert to time in ms
for i = 1:length(eegTimeStamp) 
    st = (i-1)*Nblocks + 1;
    ed = i*Nblocks;
    tsFirst = eegTimeStamp(i);
    tsEEG(st:ed) = ts + tsFirst; %add first timestamp
end

if strcmp(fn(1:3),'pie') == 1 %zoe's data
    
    Nprobes = 4; %number of probes
    Nwires = 4; %wires per probe

    channels = nan(1,Nprobes);
    %Nchannels = length(channels);
    for groupI = 1:Nprobes
        chST = (groupI-1)*Nwires + 1; %starting and ending channel
        chED = groupI*Nwires;
        chan = chST:chED; %channels
        NartOK = zeros(1,Nwires); %amount of good data
        for chI = 1:Nwires
            ch = chan(chI);
            art = signalOK{ch};
            NartOK(chI) = sum(art)/length(art);
        end
        [~,k] = max(NartOK);
        channels(groupI) = chan(k);
    end
    eegData = eegData(channels,:);
    signalOK = signalOK(channels);
else %fraser's data

    Nprobes = 8; %number of probes
    Nwires = 4; %wires per probe

    mouseID = strsplit(fn,'_'); %M16_
    mouseID = mouseID{1}; %M16
    mouseID = mouseID(2:end); %16
    mouseID = str2double(mouseID);

    disp(mouseID)

    if mouseID >= 16
        %find channels - channels are in group of 4, max 32 channels
        %go through each group and select channel with least artifacts

        channels = nan(1,Nprobes);
        %Nchannels = length(channels);
        for groupI = 1:Nprobes
            chST = (groupI-1)*Nwires + 1; %starting and ending channel
            chED = groupI*Nwires;
            chan = chST:chED; %channels
            NartOK = zeros(1,Nwires); %amount of good data
            for chI = 1:Nwires
                ch = chan(chI);
                art = signalOK{ch};
                NartOK(chI) = sum(art)/length(art);
            end
            [~,k] = max(NartOK);
            channels(groupI) = chan(k);
        end
        eegData = eegData(channels,:);
        signalOK = signalOK(channels);
    else
        %{
        M1	P11	EEG_0
        M3	P6	EEG_0
        M9	P9	my notes say that P9 should be EEG_0
            P10	my notes say that P10 should be EEG_1, but I don?t see EEG_1 in the bpf
            P11	n/a
        M10	P5	EEG_0
        M11	P11	EEG_1
        %}

        if mouseID == 1
            channels = 1*ones(1,Nprobes);
        elseif mouseID == 3
            channels = 1*ones(1,Nprobes);
        elseif mouseID == 9
            %channels = 1*ones(1,Nprobes);
            channels = [1,1,1,1,1,2,1,1];
        elseif mouseID == 10
            channels = 1*ones(1,Nprobes);
        elseif mouseID == 11
            channels = 2*ones(1,Nprobes);
        end

        eegData = eegData(channels,:);
        signalOK = signalOK(channels);
    end
end

%convert from uint32
tetrodeChannel = double(tetrodeChannel);
tetrodeUnit = double(tetrodeUnit);
tetrodeTimestamp = double(tetrodeTimestamp)/ 10; %/10 = ms, /5 = 0.5ms = eegFS
%tetrodeTimestamp = round(tetrodeTimestamp); %round to nearest ms

%find unique pairs of channel and unit
ChannelUnit = cat(2,tetrodeChannel,tetrodeUnit);
ChannelUnit = unique(ChannelUnit,'rows'); %get index of ch in unique channels
%exclude cluster 0
k = ChannelUnit(:,2) ~= 0;
ChannelUnit = ChannelUnit(k,:);

uCH = unique(tetrodeChannel); %get unique channels

resultsPH = cell(1,size(ChannelUnit,1));
resultsTS = cell(1,size(ChannelUnit,1));

%find phases of spikes
for pI = 1:size(ChannelUnit,1)
    ch = ChannelUnit(pI,1);  %channel
    u = ChannelUnit(pI,2); %unit

    %get index of ch in unique channels
    kch = uCH == ch;

    %EEG channels
    eeg = eegData(kch,:);
    art = signalOK{kch};

    %get timestamps
    k = tetrodeChannel == ch & tetrodeUnit == u;

    %timestamps of spikes (in EEG samples)
    ts = tetrodeTimestamp(k);

    %find indexes of tetrode timestamps in EEG timestamps
    edges = [-Inf, mean([tsEEG(2:end); tsEEG(1:end-1)]), +Inf];
    tsx = discretize(ts, edges); %sample indexes of EEG

    %noise free
    kOK = art(tsx);
    tsx = tsx(logical(kOK));
    if isempty(tsx); continue; end;

    phases = cell(Nbands,1);

    %filter in bands of interest
    parfor bandI = 1:Nbands

        psi = psiBands{bandI};

        %convolution of wavelet and signal
        c = conv(eeg,psi);
        %fix start and end
        N = round((length(psi)-1)/2);
        c = c(N:length(c)-N); 
        if length(c) > size(eeg,2)
            c = c(1:size(eeg,2));
        end

        phase = angle(c);

        ph = phase(tsx);

        phases{bandI} = ph;

    end

    phases = cell2mat(phases);

    resultsPH{pI} = phases;
    resultsTS{pI} = tsx;
end

save([drOUT fn],'resultsPH','resultsTS','ChannelUnit','bands');

