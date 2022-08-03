function plotMCSconnectivity(sttcs)
addpath('/home/phornauer/Git/phornauer/connectionmatlab','/home/phornauer/Git/GM-MEA/utils/2019_03_03_BCT', '/home/phornauer/Git/phornauer/Dualmode')

for f = 1:length(sttcs) % number of files
   figure('Color','w');
   tiledlayout('flow')
   for w = 1:length(sttcs{f})
       nexttile
       communityPlot(sttcs{f}{w},0)
   end
end