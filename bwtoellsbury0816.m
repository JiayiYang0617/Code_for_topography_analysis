clear
close all
%% read lidar DEM for blackwater
HAT= 3.2; 
MSL = 0.2; %for all sites

%Select a lidar tile
[file,path,indx] = uigetfile('*.tif');
if isequal(file,0)
   disp('User selected Cancel')
else
   filepath = fullfile(path, file);
   disp(['User selected ', filepath,... 
         ' and filter index: ', num2str(indx)])
end

%Display metadata
info = georasterinfo(filepath);
disp(info)
[Z1,R] = readgeoraster(filepath,'OutputType','double');

%Setup spatial coordinates for DEM
cellsize = R.CellExtentInWorldX;  %assume cells are square so just use x
X = R.XWorldLimits(1):cellsize:R.XWorldLimits(2)-1;
Y = R.YWorldLimits(1):cellsize:R.YWorldLimits(2)-1;

%Trim elevation range to threshold high ground and subtidal areas
Z1(Z1 > HAT) = HAT;
Z1(Z1 < MSL) = MSL;

%LOAD THE MASK
load('Polygon2toellsbury.mat');
Z2=Z1;
Z2(bw1>0)=nan;

%Visualise DEM raster as georeferenced image
[Xg,Yg] = meshgrid(X,Y);
figure(1)
imagesc(X,Y,Z2)
colormap gray %change color to gray
axis equal

%normalise image [0 - 1 scale]
%Z3 = Z2 - min(min(Z2));
%Z3 = Z3 ./ max(max(Z3));

Z3 = gray2ind(Z2,100);
figure(2)
imagesc(X,Y,Z3)

cmap = colormap;
Z4 = ind2rgb(Z3,cmap);
figure(3)
imagesc(X,Y,Z4)

xticks([])
yticks([])


%% SEGMENTATION
% Define thresholds for Salt pans
sliderBW1 = (Z4(:,:,1) <= 0.619) & (Z4(:,:,2) >= 0.200)...
    &(Z4(:,:,3) <= 0.860);
% Define thresholds for Saltmarsh
sliderBW2 = (Z4(:,:,1) >= 0.619);
% Define thresholds for channel
sliderBW3 = (Z4(:,:,3) >= 0.860);

% Initialize output masked image based on input image.
Z5 = Z4;
Z6 = Z4;
Z7 = Z4;

% Set background pixels where BW is false to zero.
Z5(repmat(~sliderBW1,[1 1 3])) = 0;
Z6(repmat(~sliderBW2,[1 1 3])) = 0;
Z7(repmat(~sliderBW3,[1 1 3])) = 0;

figure
imshow(Z5);
figure
imshow(Z6);
figure
imshow(Z7)
%% Absolute site area
site=Z2;
site(isnan(site))=0;
site(site>0.1)=1;
siteProps= regionprops("table",site,"Area");
sitearea=sum(siteProps.Area);
AbsoluteAreaHa=sprintf('% 2.2f%',sitearea*0.0001);
%% Pans:area and mean elevation
PansColor = rgb2ind(Z5,cmap);
Pans=Z1;
Pans(PansColor<0.2)=0;
Pans(PansColor>=0.2)=1;
PansProps = regionprops("table",Pans,"Area");
PansArea = sum(PansProps.Area);
PansArea=(PansArea/sitearea);
PansAR=sprintf('%2.1f%%', PansArea*100);

MeanPansDEM=mean(mean(Z1(PansColor>0)));
PansElevation=sprintf('%2.2f%', MeanPansDEM);

%% Saltmarsh:area and mean elevation
SaltmarshColor = rgb2ind(Z6,cmap);
Saltmarsh=Z1;
Saltmarsh(SaltmarshColor<0.2)=0;
Saltmarsh(SaltmarshColor>=0.2)=1;
SaltmarshProps = regionprops("table",Saltmarsh,"Area");
SaltmarshArea = sum(SaltmarshProps.Area);
SaltmarshArea=(SaltmarshArea/sitearea);
saltmarshAR=sprintf('%2.1f%%', SaltmarshArea*100);

MeanSaltmarshDEM=mean(mean(Z1(SaltmarshColor>0)));
SaltmarshElevation=sprintf('%2.2f%', MeanSaltmarshDEM);
%% Channel: area
ChannelColor = rgb2ind(Z7,cmap);
Channel=Z1;
Channel(ChannelColor<0.2)=0;
Channel(ChannelColor>=0.2)=1;
ChannelProps = regionprops("table",Channel,"Area");
ChannelArea = sum(ChannelProps.Area);
ChannelArea=(ChannelArea/sitearea);
ChannelAR=sprintf('%2.1f%%', ChannelArea*100);
%% Define ROI

%Interactively define a polygon (digitise, right-click > create mask)
I = mat2gray(Z1);  %we need a grayscale image to use ROIPOLY function
[BW, xi,yi] = roipoly(I);

% *** Note - would be good idea to add code to save and then re-use
% polygons to avoid interactively creating them every time the script is
% run. Maybe create once, save, and then used saved one. Code could be made
% to check for existing polygon.
%Crop image using ROI to reduce memory use
bw1=imcomplement(BW);
save('Polygon2toellsbury.mat','xi','yi','bw1');
load('Polygon2toellsbury.mat');

%% Save all Zdata
Zdata_bwtoellsbury.DEM=Z1;
Zdata_bwtoellsbury.maskedDEM=Z2;
Zdata_bwtoellsbury.ind=Z3;
Zdata_bwtoellsbury.rgb=Z4;
Zdata_bwtoellsbury.pans=Z5;
Zdata_bwtoellsbury.saltmarsh=Z6;
Zdata_bwtoellsbury.channel=Z7;
Zdata_bwtoellsbury.SiteArea=AbsoluteAreaHa;
Zdata_bwtoellsbury.pansarea=PansAR;
Zdata_bwtoellsbury.meanpansdem=PansElevation;
Zdata_bwtoellsbury.saltmarsharea=saltmarshAR;
Zdata_bwtoellsbury.meansaltmarshdem=SaltmarshElevation;
Zdata_bwtoellsbury.channelarea=ChannelAR;
save('Zdata_bwtoellsbury.mat','-v7.3','Zdata_bwtoellsbury')