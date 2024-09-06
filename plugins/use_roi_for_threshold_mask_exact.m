function a_use_roi_threshold_and_open(DATA,ROI,AddInfo)
% This Aedes plugin creates invert MASK based on ROI data. Uses first ROI.

% This function is written for Aedes
%
% (c) 07/2012, Mikko Nissi, nissi@cmrr.umn.edu
%
%
% This program may be used under the terms of the GNU General Public
% License version 2.0 as published by the Free Software Foundation
% and appearing in the file LICENSE.TXT included in the packaging of
% this program.
%
% This program is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
% WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

% grab data (use only current volume)
if length(ROI)==0
    disp('No ROI.');
    return
end

emelent3(:,:,1)=[...
    0 0 0;...
    0 1 0;...
    0 0 0];

emelent3(:,:,2)=[...
    0 1 0;...
    1 1 1;...
    0 1 0];

emelent3(:,:,3)=[...
    0 0 0;...
    0 1 0;...
    0 0 0];


emelent(:,:,1)=[...
     0     0     0     0     0;...
     0     1     1     1     0;...
     0     1     1     1     0;...
     0     1     1     1     0;...
     0     0     0     0     0];
 
emelent(:,:,2)=[...
     0     0     1     0     0;...
     0     1     1     1     0;...
     1     1     1     1     1;...
     0     1     1     1     0;...
     0     0     1     0     0];

emelent(:,:,3)=[...
     0     1     1     1     0 ;...
     1     1     1     1     1 ;...
     1     1     1     1     1 ;...
     1     1     1     1     1 ;...
     0     1     1     1     0 ];

emelent(:,:,4)=[...
     0     0     1     0     0;...
     0     1     1     1     0;...
     1     1     1     1     1;...
     0     1     1     1     0;...
     0     0     1     0     0];

emelent(:,:,5)=[...
     0     0     0     0     0;...
     0     1     1     1     0;...
     0     1     1     1     0;...
     0     1     1     1     0;...
     0     0     0     0     0];
 
emelent=logical(emelent);
 
% St.Dev. threshold factor
%sdfac=2;
sdfac=1;

% grab and smooth data
dat=DATA{1}.FTDATA;
%dat=smooth3(dat,'box',3);

% grab values from ROI
roi=dat(ROI.voxels{1});

% no forgiving, use ROI avg for threshold
datlim=mean(roi); %-sdfac*std(roi);



% threshold
dat(dat<datlim)=0;
dat(dat>=datlim)=1;

% stuff this mask back to ROI variable to be saved later..
ROI.voxels{1}=logical(dat);

% 
% 
% % do some closing operations

%st=logical(dat);

% Close once: erode & dilate
%st=dilate3d(st,emelent3);
%st=erode3d(st,emelent3);
%st=erode3d(st,emelent3);

%ROI.voxels{1}=st;

% st=dilate3d(st,emelent);
% st=erode3d(st,emelent);
% st=dilate3d(st,emelent);
% st=erode3d(st,emelent);
% %st=dilate3d(st);
% %st=erode3d(st);
% %st=dilate3d(st);
% %st=erode3d(st);
% %st=dilate3d(st);
% %st=erode3d(st);
% 
% % invert
% st=logical(1-st);
% % and dilate..
% st=dilate3d(st,emelent);
% st=dilate3d(st,emelent);
% st=dilate3d(st,emelent);
% 

% check if name was given, then don't ask
if (nargin<4) || isempty(savefname)
% Prompt for file name
[fname,fpath,findex]=uiputfile({'*.mat;*.Mat;*.roi;*.ROI;*.Roi*.MAT;*.*',...
                    'roi-Files (*.mat;*.Mat;*.MAT;*.roi;*.ROI;*.Roi)';...
                    '*.*','All Files (*.*)'},...
                               'Save roi-file',[DATA{1}.HDR.fpath, ...
                    'mask_ROI.roi']);
if isequal(fname,0) || isequal(fpath,0)
  return
end
else
    % save-filename was given, make compatible with the rest of the
    % function
    [fpath,fname,fext]=fileparts(savefname);
    fname=[fname fext];
end

% Calculate the map
[fp,fn,fe]=fileparts([fpath,fname]);


save(fullfile(fpath,fname),'-mat','ROI');


    
