% clear all

% paths
addpath(genpath('/home/phornauer/Git/STTC'));
addpath(genpath('/home/phornauer/Git/MCS'));
addpath(genpath('/home/phornauer/Git/GM-MEA/utils/NCCToolboxV1'));

base_dir = '/home/phornauer/Data/Clemens/Recordings';
save_folder = '/home/phornauer/Data/Clemens/Connectivity';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Infer STTC w probabilistic threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% parameter for sttc 
params.fs = 1000000; % sampling rate
params.max_lag = 0.05; % in seconds, the delay time to look for conincidental spikes
params.iter = 100; % iterations used for surrogate generation
params.surr_v = 'jitter'; % method for surrogate generation
params.pct_thresh = 95; % percentile threshold
params.min_rate = 0;%0.01; % Hz
params.max_rate = 1000;%10; % Hz
params.rec_length = 12000; % in seconds

% find the spike data
all_files = dir(fullfile(base_dir,'*.xlsx'));
sttcs = cell(1,length(all_files));

for i = 1:length(all_files)
%     try
        file_sttc = cell(1,24);
        file_name = fullfile(all_files(i).folder, all_files(i).name);
        file_id = strsplit(all_files(i).name, '.');
        save_path = fullfile(save_folder, file_id{1});
        fprintf(['Processing: ' file_name '\n']);
        [spikes, ~, ~] = spiketrainesFromMCSExcel(file_name,params.fs);
        for w = 1:length(spikes)
            fprintf('Well %i\n',w)
            if ~exist(save_path,'dir')
                mkdir(save_path)
            end
            save_file_name = fullfile(save_path,sprintf('well_%i',w));
            [sttc_results] = sttc_analysis(spikes{w},save_file_name,params);
            file_sttc{w} = sttc_results.sttc_bu;
        end
        sttcs{i} = file_sttc;
%     catch ME
%         fprintf('%s %s \n', ME.identifier, ME.message)
%     end
end
