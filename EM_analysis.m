function EM_analysis
%% Adam Tyson | 29/03/2018 | adam.tyson@icr.ac.uk
% loads EM image, allows manual segmentation of >=1 cell, and then 
% performs automatic analysis

%% RETURNS
% Area - number of pixels in segmented image
% MajorAxisLength & MinorAxisLength - of normalised ellipse
% Convex area - convex bounding shape area
% Perimeter - of segmented shape
% AxisRatio - ratio of Minor to Major axes
% concave ratio - pixels in segmented image divided by convex area 
%               - measure of protrusions etc
% PerimAreaRatio - another measure of protrusions
% Mean - mean intensity of cell (probably not much use for EM)
% Entropy - measure of texture (how different pixels are)

%% TO DO
% add option to only analyse certain images
% add more texture measurements

vars=getVars;

cd(vars.directory) 
files=dir('*.tif'); % all tif's in this folder
numImages=length(files);
imCount=0;

for file=files' % go through all images
     imCount=imCount+1; 
%% Load images and separate cells
im.raw{imCount}= imread(file.name);
[im.binary{imCount}, cellNum(imCount)] = manSeg(im.raw{imCount});
im.cell_indiv{imCount} = maskCell(im.raw{imCount},im.binary{imCount}, cellNum(imCount));
end

%% analyse cells
tic
progressbar('Analysing images') % Init prog bar
count=0;
imCount=0;
for iInd=1:numImages % go through all images
    count=count+1; 

     for cInd=1:length(im.cell_indiv{iInd})    
         imCount=imCount+1;
         cellName{imCount}=[files(iInd).name '__cell_' num2str(cInd)];
         
         paramTmp = EM_describe(im.cell_indiv{1,iInd}{cInd,1});
         cellParams(imCount,:)=cell2mat(struct2cell(paramTmp));
     end 
     
     % progress bar
     frac1 =count/numImages;
     progressbar(frac1)
end

save_res(cellParams, cellName, vars)
        
disp(['Time elapsed: ' num2str(toc) ' seconds'])

end

%% Internal functions
function [binaryImages, cellNum]=manSeg(image)

    scrsz = get(0,'ScreenSize');
    imSize=size(image);
    dispScale=(scrsz(4)/imSize(1))*0.8;
    screenSize=[10 10 dispScale*imSize(2) dispScale*imSize(1)];

    figure('position', screenSize,'Name','Manually segment cells')
    imagesc(image)
    colormap gray

    continueSeg=1;
    cellNum=1;
    
    while continueSeg==1
        hFH = imfreehand(); % manually segment
        tmpBin = hFH.createMask(); % make binary image
        repSeg = questdlg('Redo last segmentation?',...
            'Error catch','Yes','No','No'); % yes, no and default 
        
        if strcmp(repSeg, 'No')  
            binaryImages(:,:,cellNum) = tmpBin;
            finSeg = questdlg('All cells segmented?',...
                'Next image','Yes','No','No');

            if strcmp(finSeg, 'No') 
                cellNum = cellNum+1;
            elseif strcmp(finSeg, 'Yes')
                continueSeg = 0;
                close all
            end
            
        end
    end
end

function imageCrop=deleteZeros(image)
    image_max = max(image, [], 3);
        for z=1:size(image,3)
            imagetmp=image(:,:,z);
            imagetmp( all(~image_max,2), :) = []; % remove zero rows
            imagetmp( :, all(~image_max,1)) = []; % remove zero columns
            imageCrop(:,:,z)=imagetmp;
        end
end

function cell_indiv = maskCell(rawimage, binaryImages, cellNum)
    cell_indiv = cell(cellNum, 1) ;
   
    for cellNum=1:cellNum
        C0_indv_tmp=rawimage.*uint8(binaryImages(:,:,cellNum));       
        imageCrop_C0=deleteZeros(C0_indv_tmp);
        cell_indiv{cellNum}=imageCrop_C0;
    end
  
end

function vars=getVars
    vars.directory = uigetdir('', 'Choose directory containing images');
    
    vars.stamp=num2str(fix(clock)); % date and time 
    vars.stamp(vars.stamp==' ') = '';%remove spaces
end

function save_res(cellParams, cellName, vars)
    fill_side=zeros(size(cellParams,1), 1);
    cellParams=[fill_side cellParams];
    cellParams_cell=num2cell(cellParams);
    
    for i=1:length(cellName)
        cellParams_cell{i,1}=cellName{i};
    end
    
    tableNames{1,1}='Cell_Parameter';
    tableNames{1,2}='Area';
    tableNames{1,3}='Major_axis_length';
    tableNames{1,4}='Minor_axis_length';
    tableNames{1,5}='Convex_area';
    tableNames{1,6}='Perimeter';
    tableNames{1,7}='Minor_major_axis_ratio';
    tableNames{1,8}='Area_concave_area_ratio';
    tableNames{1,9}='Perimeter_area_ratio';
    tableNames{1,10}='Mean_intensity';
    tableNames{1,11}='Entropy';

    results_Table=cell2table(cellParams_cell);
    results_Table.Properties.VariableNames = tableNames;
    writetable(results_Table, ['Results_' vars.stamp '.csv'])
end
