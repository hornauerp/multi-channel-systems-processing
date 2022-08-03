function [spikes, well_labels, channel_labels] = spiketrainesFromMCSExcel(file_path,sampling_rate)
data = readtable(file_path);
well_labels = unique(data.WellLabel);
channel_labels = unique(data.ChannelLabel);
spikes = cell(1,24);
for w = 1:length(well_labels)
    well_data = data(ismember(data.WellLabel,well_labels{w}),:); %Select rows with corresponding well label
    well_spikes = cell(1,12);
    for c = 1:length(channel_labels)
        well_spikes{c} = well_data.Timestamp(ismember(well_data.ChannelLabel,channel_labels{c}))/sampling_rate;
    end
    spikes{w} = well_spikes;
end