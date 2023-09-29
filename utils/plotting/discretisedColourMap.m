function rgbs = discretisedColourMap(colourmap,flip,numBins)
% colourmap='winter'; flip=true; numBins=3;

if flip
    temp0 = flipud(eval(colourmap));
else
    temp0 = eval(colourmap);
end
[~,temp1] = discretize([1,numBins],size(temp0,1));
temp2 = discretize(1:numBins,temp1);
rgbs = temp0(temp2,:);