addpath(genpath('/home/phornauer/Git/phornauer/connectionmatlab/FMAToolbox/Helpers'))
%%
ccg_params.binSize = .001; % .5ms 
ccg_params.duration = .1; %50ms, corresponds to the length of one side of ccg
ccg_params.epoch = [0 inf]; %whole session
ccg_params.conv_w = .010/ccg_params.binSize;  % 10ms window (gaussian convolution)
ccg_params.alpha = 0.001; %high frequency cut off, must be .001 for causal p-value matrix -> that's the line from the original script but you can use whatever... :)
ccg_params.Fs = 1/1000000;
% find the spike data
all_files = dir(fullfile(base_dir,'*.xlsx'));
sig_cons = [];
for i = 1:length(all_files)
%     try
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
            [sig_con, ccg_vec, Bounds, ccg_vec_inh, sig_con_inh,ccgR, Pval] = ccgMCS(spikes{w},ccg_params);
            if ~isempty(sig_con)
                sig_cons = vertcat(sig_cons,[i w sig_con]);
            end
        end
        
%     catch ME
%         fprintf('%s %s \n', ME.identifier, ME.message)
%     end
end