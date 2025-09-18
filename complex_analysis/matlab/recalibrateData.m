function [ptEstimateNew,errorBarNew,combinedPtEstIntNew,...
    combinedSensIntNew] = recalibrateData(ptEstimate,errorBar,...
    combinedPtEstInt,combinedSensInt,r_1,r_2,t_1,t_2,deltaF)

combinedSensIntNew = combinedSensInt./(r_1.*r_1.*r_2.*r_2);

l4l3 = sum(deltaF.*combinedSensInt)./...
    sum(deltaF.*combinedSensInt./(r_1.*r_1.*r_2.*r_2));

combinedPtEstIntNew =  exp(1i.*(t_2-t_1)) .* l4l3 .* ...
        combinedPtEstInt ./(r_1.*r_2);
    
ptEstimateNew=sum(deltaF.*2.*real(combinedPtEstIntNew));
errorBarNew = sqrt(1/sum(combinedSensIntNew.*deltaF));

return
