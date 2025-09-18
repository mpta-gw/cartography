% combineResultsFromMultiplePairs - script to combine the final post-processed 
% results from multiple detector pairs. 
%
% Input:
%	outputDir - directory where all post processing results are located
%		    and to which all the combined results will be written.
%	ifos      - a string containing all the names of the ifos that you
%		    wish to combine, of the form 'H1L1V1', i.e. with no
%		    separator. The script will form this into pairs; if a
%		    pair does not exist, it will skip it (with a warning).
%
% Written by Emma Robinson
% Contact elr@star.sr.bham.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function combineResultsFromMultiplePairs(outputDir,ifos)

% ifos is a list of all detectors, with no separator, eg. H1V1L1, so must have
% a length which is even

if mod(length(ifos),2)
	error('All IFOs should have a two character name')
end

if ~exist(outputDir,'dir')
	error('Output directory does not exist')
end

numIFOS = length(ifos)/2;

numPairs = numIFOS * (numIFOS-1) / 2;

pairnum=1;

pair=[];
ptEstimate = zeros(1,numPairs);
errorBar = zeros(1,numPairs);
existent=[];
for ii=1:numIFOS-1
	for jj=ii+1:numIFOS
		pair = [pair;[ ifos((2*ii)+[-1:0]) , ifos((2*jj)+[-1:0]) ]];

		ptEstIntFile = [outputDir , '/' , pair(pairnum,:) , '_ptEstIntegrand.dat'];
		sensIntFile  = [outputDir , '/' , pair(pairnum,:) , '_sensIntegrand.dat' ];

		if exist(ptEstIntFile,'file') && exist(sensIntFile,'file')
			if ii==1 && jj==2 %first pass
				[temp, flow, deltaF, errorBar(pairnum)] = readCombinedSensIntFromFile(sensIntFile);
				sensInts = zeros(length(temp),numPairs);
				sensInts(:,pairnum) = temp;
				[temp, flow, deltaF, errorBar(pairnum)] = readCombinedPtEstIntFromFile(ptEstIntFile);
                                ptEstInts = zeros(length(temp),numPairs);
                                ptEstInts(:,pairnum) = temp;

			else
				[sensInts(:,pairnum), flow, deltaF, errorBar(pairnum)] = readCombinedSensIntFromFile(sensIntFile);
				[ptEstInts(:,pairnum), flow, deltaF, errorBar(pairnum)] = readCombinedPtEstIntFromFile(ptEstIntFile);
			end

			ptEstimate(pairnum) = deltaF*2*sum(real(ptEstInts(:,pairnum)));
			
			existent=[existent,pairnum];
		else
			pair(pairnum,:) = [ ifos((2*jj)+[-1:0]) , ifos((2*ii)+[-1:0]) ];

			ptEstIntFile = [outputDir , '/' , pair(pairnum,:) , '_ptEstIntegrand.dat'];
	                sensIntFile  = [outputDir , '/' , pair(pairnum,:) , '_sensIntegrand.dat' ];

        	        if exist(ptEstIntFile,'file') && exist(sensIntFile,'file')
                	        if ii==1 && jj==2 %first pass
                        	        [temp, flow, deltaF, errorBar(pairnum)] = readCombinedSensIntFromFile(sensIntFile);
                                	sensInts = zeros(length(temp),numPairs);
	                                sensInts(:,pairnum) = temp;
        	                        [temp, flow, deltaF, errorBar(pairnum)] = readCombinedPtEstIntFromFile(ptEstIntFile);
                	                ptEstInts = zeros(length(temp),numPairs);
                        	        ptEstInts(:,pairnum) = temp;

	                        else
        	                        [sensInts(:,pairnum), flow, deltaF, errorBar(pairnum)] = readCombinedSensIntFromFile(sensIntFile);
                	                [ptEstInts(:,pairnum), flow, deltaF, errorBar(pairnum)] = readCombinedPtEstIntFromFile(ptEstIntFile);
                        	end

	                        ptEstimate(pairnum) = deltaF*2*sum(real(ptEstInts(:,pairnum)));
	                        existent=[existent,pairnum];

			else
				warning(['Data is not present for pair ' pair(pairnum,:) '!']);
			end
		end

		pairnum=pairnum+1;
	end
end
%get rid of non-existent pairs
%existent = ~isnan(1./errorBar);
errorBar=errorBar(existent);
ptEstimate=ptEstimate(existent);
ptEstInts=ptEstInts(:,existent);
sensInts=sensInts(:,existent);
pair=pair(existent,:);
nn=size(pair);
numPairs=nn(1);
clear nn
%combined results

combinedPtEstimate = ( sum((ptEstimate.*(errorBar.^(-2)))))./(sum(errorBar.^(-2)));

combinedErrorBar = 1 / sqrt( sum( (errorBar.^(-2)) ));

combinedSensInt = sum(sensInts,2);

numFreqs = length(combinedSensInt);
freqs = flow + (deltaF*[0:numFreqs-1]);

save([outputDir, '/' , 'combinedResults_' , ifos ,'.mat']);

fignum=1;

%first plot - pt est and error bar



figure(fignum)
h=errorbar([1:numPairs+1],[ptEstimate,combinedPtEstimate],[errorBar,combinedErrorBar],'.-');
set(h,'LineStyle','none')
set(h,'MarkerSize',15)
xt=[0:1:numPairs+2];
xtick=cell(1,numPairs+3);

xtick{1}='';
for ii=2:numPairs+1
        xtick{ii}=pair(ii-1,:);
end
xtick{numPairs+2}=ifos;

xlim([0 numPairs+2]);
set(gca,'XTick',xt);
set(gca,'XTickLabel',xtick);
set(gca,'XGrid','on')
set(gca,'YGrid','on')
set(gca,'FontSize',12)
ylabel('h_{100}^2\Omega_{gw}(f_R)','FontSize',14)
title('Point Estimates','FontSize',16')

print(fignum,'-depsc2',[outputDir, '/' , 'combinedPtEstErrBar_' ifos])
print(fignum,'-dpng',[outputDir, '/' , 'combinedPtEstErrBar_' ifos])
saveas(fignum,[outputDir, '/' , 'combinedPtEstErrBar_' ifos '.fig'])

close(fignum)
fignum=fignum+1;
%second figure - sensitivity integrands
figure(fignum)
fig_handle = plot( freqs ,[sensInts, combinedSensInt]);

set(gca,'FontSize',12);
xlabel('frequency (Hz)','FontSize',14);
ylabel('I(f) (Hz^2 / strain^4)','FontSize',14);
title('Stochastic Sensitivity Integrand','FontSize',14);
set(gca,'XGrid','on')
set(gca,'YGrid','on')
for ii=1:numPairs
	legcell{ii}=pair(ii,:);
end
legcell{ii+1}=ifos;
leg=legend(legcell);
set(leg,'FontSize',14);
print(fignum,'-depsc2',[outputDir, '/' , 'combinedSensitivityIntegrand_' ifos])
print(fignum,'-dpng',[outputDir, '/' , 'combinedSensitivityIntegrand_' ifos])
saveas(fignum,[outputDir, '/' , 'combinedSensitivityIntegrand_' ifos '.fig'])
close(fignum)
