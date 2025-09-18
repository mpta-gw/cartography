load('S5')
[d,RA,DECL,dOmg,U]=diagPixel(fisher_opt,res);
x_radio=U*x_opt;
p_radio=x_radio./d;
save('S5_radio.mat','d','RA','DECL','dOmg','U','x_radio','p_radio');

plotMapAitoff(p_radio,360,181,-1);
print('-depsc2','radio_S5.eps')
print('-djpeg','radio_S5.jpeg')
