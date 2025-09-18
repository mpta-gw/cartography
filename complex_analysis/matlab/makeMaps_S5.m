load('S5.mat');

plotMapAitoff(map./sigmaPix,360,181,-1);
print('-depsc2','SNR_S5.eps')
print('-djpeg','SNR_S5.jpeg')

plotMapAitoff(sigmaPix,360,181,-1);
print('-depsc2','sigma_S5.eps')
print('-djpeg','sigma_S5.jpeg')

plotMapAitoff(ul_map',360,181,-1);
print('-depsc2','ul_S5.eps')
print('-djpeg','ul_S5.jpeg')

plotPlmMap(p_opt);
print('-depsc2','clean_S5.eps')
print('-djpeg','clean_S5.jpeg')

plotPlmMap(x_opt);
print('-depsc2','dirty_S5.eps')
print('-djpeg','dirty_S5.jpeg')

