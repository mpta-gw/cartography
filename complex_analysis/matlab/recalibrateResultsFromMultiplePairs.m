% recalibrateResultsFromMultipleDetectors - script to recalibrate the
% results from multiple detectors
%
% Input:
%   originalMatFile        - .mat file in which original results were saved
%   outputDir              - directory into which .mat file containing
%                            new results will be saved
%   recalibrationFileNames - cell array of names of files in which the 
%                            functions to convert between versions are 
%                            written. If a detector does not need to be
%                            recalibrated, the file name should be []
function recalibrateResultsFromMultiplePairs(originalMatFile,...
            outputDir,recalibrationFileNames)
        
% load original data
originalData=load(originalMatFile);

% read the recalibration data
for ii=1:originalData.numIFOS
    if isempty(recalibrationFileNames{ii})
        r{ii}=ones(size(originalData.ptEstInts,1),1);
        t{ii}=zeros(size(originalData.ptEstInts,1),1);
    else
        tmp=load(recalibrationFileNames{ii});
        [fc,i1,i2]=find_freqcut(originalData.freqs(1),originalData.freqs(end)+originalData.deltaF,tmp(:,1));
        r{ii} = interp1(fc,tmp(i1:i2,2),originalData.freqs).';
        t{ii} = interp1(fc,tmp(i1:i2,3),originalData.freqs).';
    end
end

%initialise

ptEstimateNew=zeros(size(originalData.ptEstimate));
errorBarNew=zeros(size(originalData.errorBar));
ptEstIntsNew=zeros(size(originalData.ptEstInts));
sensIntsNew=zeros(size(originalData.sensInts));

kk=1;
for ii=1:originalData.numIFOS-1
    for jj=ii+1:originalData.numIFOS
        
        if originalData.pair(kk,:)==[ originalData.ifos((2*ii)+[-1:0]) , ...
                originalData.ifos((2*jj)+[-1:0]) ]
            fprintf('pair = %s\n',originalData.pair(kk,:));

            
            [ptEstimateNew(kk),errorBarNew(kk),ptEstIntsNew(:,kk),...
                sensIntsNew(:,kk)] = ...
                recalibrateData(originalData.ptEstimate(kk),...
                originalData.errorBar(kk),originalData.ptEstInts(:,kk),...
                originalData.sensInts(:,kk),r{ii},r{jj},...
                t{ii},t{jj},originalData.deltaF);            
            
            fprintf('original:\t%f+/-%f\n',originalData.ptEstimate(kk),...
                originalData.errorBar(kk));
            fprintf('new:\t%f+/-%f\n\n',ptEstimateNew(kk),errorBarNew(kk));

            kk=kk+1;
        else
            fprintf('missing\n');
        end
    end
end

%combined results

combinedPtEstimateNew = ( sum((ptEstimateNew.*(errorBarNew.^(-2)))))./(sum(errorBarNew.^(-2)));


combinedPtEstIntNew = zeros(size(ptEstIntsNew,1),1);
for ii=1:originalData.numPairs
    combinedPtEstIntNew = combinedPtEstIntNew + ...
        (((ptEstIntsNew(:,ii)./(errorBarNew(ii).^2))/(errorBarNew(ii).^-2)));
end
if ~isfield(originalData,'combinedPtEstInt')
    originalData.combinedPtEstInt= zeros(size(originalData.ptEstInts,1),1);
    for ii=1:originalData.numPairs
        originalData.combinedPtEstInt= originalData.combinedPtEstInt+ ...
            (((originalData.ptEstInts(:,ii)./...
            (originalData.errorBar(ii).^2))/...
            (originalData.errorBar(ii).^-2)));
    end
end


combinedErrorBarNew = 1 / sqrt( sum( (errorBarNew.^(-2)) ));

combinedSensIntNew = sum(sensIntsNew,2);

save([outputDir, '/' , 'recalibratedResults_' , originalData.ifos ,'.mat']);

fprintf('COMBINED:\n')
fprintf('original:\t%f+/-%f\n',originalData.combinedPtEstimate,...
    originalData.combinedErrorBar);
fprintf('new:\t%f+/-%f\n\n',combinedPtEstimateNew,combinedErrorBarNew);

%first plot - pt est and error bar

fignum=1;

figure(fignum)
h=errorbar([1:originalData.numPairs+1],[originalData.ptEstimate,originalData.combinedPtEstimate],...
    [originalData.errorBar,originalData.combinedErrorBar],'.-');
set(h,'LineStyle','none')
set(h,'MarkerSize',15)
hold on
h2=errorbar([1:originalData.numPairs+1],[ptEstimateNew,combinedPtEstimateNew],...
    [errorBarNew,combinedErrorBarNew],'r.-');
set(h2,'LineStyle','none')
set(h2,'MarkerSize',15)
xt=[0:1:originalData.numPairs+2];
xtick=cell(1,originalData.numPairs+3);

xtick{1}='';
for ii=2:originalData.numPairs+1
        xtick{ii}=originalData.pair(ii-1,:);
end
xtick{originalData.numPairs+2}='all';

xlim([0 originalData.numPairs+2]);
set(gca,'XTick',xt);
set(gca,'XTickLabel',xtick);
set(gca,'XGrid','on')
set(gca,'YGrid','on')
set(gca,'FontSize',12)
ylabel('h_{100}^2\Omega_{gw}(f_R)','FontSize',14)
title('Point Estimates','FontSize',16')
leg=legend('before','after');

print(fignum,'-depsc2',[outputDir, '/' , 'recalibratedPtEstErrBar_' originalData.ifos])
print(fignum,'-dpng',[outputDir, '/' , 'recalibratedPtEstErrBar_' originalData.ifos])
saveas(fignum,[outputDir, '/' , 'reacalibratedPtEstErrBar_' originalData.ifos '.fig'])

close(fignum)
fignum=fignum+1;

%second figure - sensitivity integrands
figure(fignum)
fig_handle = plot( originalData.freqs ,[sensIntsNew, combinedSensIntNew]);

set(gca,'FontSize',12);
xlabel('frequency (Hz)','FontSize',14);
ylabel('I(f) (Hz^2 / strain^4)','FontSize',14);
title('Stochastic Sensitivity Integrand','FontSize',14);
set(gca,'XGrid','on')
set(gca,'YGrid','on')
for ii=1:originalData.numPairs
	legcell{ii}=originalData.pair(ii,:);
end
legcell{ii+1}=originalData.ifos;
leg=legend(legcell);
set(leg,'FontSize',14);
print(fignum,'-depsc2',[outputDir, '/' , 'recalibratedSensitivityIntegrand_' originalData.ifos])
print(fignum,'-dpng',[outputDir, '/' , 'recalibratedSensitivityIntegrand_' originalData.ifos])
saveas(fignum,[outputDir, '/' , 'recalibratedSensitivityIntegrand_' originalData.ifos '.fig'])
close(fignum)
fignum=fignum+1;
%third figure - compare combined sensitivity integrands
figure(fignum)
fig_handle = plot( originalData.freqs ,[originalData.combinedSensInt, combinedSensIntNew]);
set(gca,'FontSize',12);
xlabel('frequency (Hz)','FontSize',14);
ylabel('I(f) (Hz^2 / strain^4)','FontSize',14);
title('Stochastic Sensitivity Integrand','FontSize',14);
set(gca,'XGrid','on')
set(gca,'YGrid','on')
for ii=1:originalData.numPairs
	legcell{ii}=originalData.pair(ii,:);
end
leg=legend('before','after');
set(leg,'FontSize',14);
print(fignum,'-depsc2',[outputDir, '/' , 'recalibratedSensitivityIntegrandCF_' originalData.ifos])
print(fignum,'-dpng',[outputDir, '/' , 'recalibratedSensitivityIntegrandCF_' originalData.ifos])
saveas(fignum,[outputDir, '/' , 'recalibratedSensitivityIntegrandCF_' originalData.ifos '.fig'])
close(fignum)

