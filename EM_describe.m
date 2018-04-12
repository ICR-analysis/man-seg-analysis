function params = EM_describe(image)
%% Adam Tyson | 29/03/2018 | adam.tyson@icr.ac.uk
% image=im.cell_indiv{1,iInd}{cInd,1};

binIm=zeros(size(image));
binIm(image>0)=1;

params = regionprops(binIm,'Area','MajorAxisLength','MinorAxisLength',...
    'ConvexArea', 'Perimeter');

params.AxisRatio=params.MinorAxisLength/params.MajorAxisLength;
params.ConcaveRatio=params.Area/params.ConvexArea;
params.PerimAreaRatio=params.Perimeter/params.Area;

params.Mean=mean(nonzeros(image(:)));
params.Entropy=entropy(nonzeros(image(:)));

end