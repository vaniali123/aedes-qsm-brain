function Out = aedes_rotateflip(DATA,squareImages)
% AEDES_ROTATEFLIP - GUI for rotating and flipping files
%   
%
% Synopsis: 
%
% Description:
%
% Examples:
%
% See also:
%

% This function is a part of Aedes - A graphical tool for analyzing 
% medical images
%
% Copyright (C) 2006 Juha-Pekka Niskanen <Juha-Pekka.Niskanen@uku.fi>
% 
% Department of Physics, Department of Neurobiology
% University of Kuopio, FINLAND
%
% This program may be used under the terms of the GNU General Public
% License version 2.0 as published by the Free Software Foundation
% and appearing in the file LICENSE.TXT included in the packaging of
% this program.
%
% This program is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
% WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.


% Public variables
H = [];
Out = [];
Dat = [];
cancel = true;

% Get file names and paths
for ii=1:length(DATA)
  Dat.fnames{ii} = DATA{ii}.HDR.fname;
  Dat.fpaths{ii} = DATA{ii}.HDR.fpath;
  Dat.fullfiles{ii} = [DATA{ii}.HDR.fpath,...
                      DATA{ii}.HDR.fname];
end
H = l_DrawGUI;

% Wait for quit
waitfor(H.FH)
if cancel
  clear H Dat cancel
  return
end

function H=l_DrawGUI()

%% Load default font and colors
%FigColor=get(0,'DefaultUicontrolBackgroundcolor');
GD=aedes_gui_defaults;
%GD.col.mainfig = FigColor;
fig_h = 305;
fig_w = 500;
fig_location = aedes_dialoglocation([fig_w,fig_h]);
fig_pos = [fig_location(1) fig_location(2) fig_w fig_h];

%% The main figure
H.FH = figure('position',fig_pos,...
              'Units','Pixel', ...
              'Name','Rotate/Flip images', ...
              'Numbertitle','off', ...
              'Tag','im_rotate_gui', ...
              'Color',GD.col.mainfig, ...
              'Toolbar','none', ...
              'Menubar','none', ...
              'DoubleBuffer','on', ...
              'DockControls','off',...
              'renderer','painters',...
              'KeyPressFcn','',...
              'resize','off',...
              'windowstyle','modal');

% Options uipanel
H.OPTUIPANEL = uipanel('parent',H.FH,...
                       'units','pixel',...
					   'backgroundcolor',GD.col.mainfig, ...
                       'position',[5 40 fig_w-10 fig_h-45]);

% Rotation uipanel
H.ROTUIGRP = uibuttongroup('parent',H.OPTUIPANEL,...
                           'units','pixel',...
                           'position',[10 135 110 105],...
                           'title','Rotation (CCW)',...
						   'backgroundcolor',GD.col.mainfig, ...
                           'SelectionChangeFcn',@l_SelectionChanged);

% Rotation Radio buttons
H.ROTNONE=uicontrol('parent',H.ROTUIGRP,...
                    'units','pixel',...
                    'style','radio',...
                    'position',[20 70 80 15],...
					'backgroundcolor',GD.col.mainfig, ...
                    'string','none',...
                    'value',1);
tmp=get(H.ROTNONE,'position');
H.ROT90=uicontrol('parent',H.ROTUIGRP,...
                  'units','pixel',...
                  'style','radio',...
                  'position',[tmp(1) tmp(2)-20 tmp(3) tmp(4)],...
                  'string',['90',char(186)],...
				  'backgroundcolor',GD.col.mainfig, ...
                  'value',0);
tmp=get(H.ROT90,'position');
H.ROT180=uicontrol('parent',H.ROTUIGRP,...
                   'units','pixel',...
                   'style','radio',...
                   'position',[tmp(1) tmp(2)-20 tmp(3) tmp(4)],...
                   'string',['180',char(186)],...
				   'backgroundcolor',GD.col.mainfig, ...
                   'value',0);
tmp=get(H.ROT180,'position');
H.ROT270=uicontrol('parent',H.ROTUIGRP,...
                   'units','pixel',...
                   'style','radio',...
                   'position',[tmp(1) tmp(2)-20 tmp(3) tmp(4)],...
                   'string',['270',char(186)],...
				   'backgroundcolor',GD.col.mainfig, ...
                   'value',0);
if ~squareImages
  set([H.ROT90,H.ROT270],'enable','off')
end

% Flip uipanel
tmp = get(H.ROTUIGRP,'position');
H.FLIPUIGRP = uibuttongroup('parent',H.OPTUIPANEL,...
                            'units','pixel',...
                            'position',[tmp(1) tmp(2)-85-15 tmp(3) 85],...
                            'title','Flipping',...
							'backgroundcolor',GD.col.mainfig, ...
                            'SelectionChangeFcn',@l_SelectionChanged);

% Flipping radio buttons
H.FLIPNONE=uicontrol('parent',H.FLIPUIGRP,...
                     'units','pixel',...
                     'style','radio',...
                     'position',[20 50 80 15],...
					 'backgroundcolor',GD.col.mainfig, ...
                     'string','none',...
                     'value',1);
tmp=get(H.FLIPNONE,'position');
H.FLIPUD=uicontrol('parent',H.FLIPUIGRP,...
                   'units','pixel',...
                   'style','radio',...
                   'position',[tmp(1) tmp(2)-20 tmp(3) tmp(4)],...
                   'string','Up/Down',...
				   'backgroundcolor',GD.col.mainfig, ...
                   'value',0);
tmp=get(H.FLIPUD,'position');
H.FLIPLR=uicontrol('parent',H.FLIPUIGRP,...
                   'units','pixel',...
                   'style','radio',...
                   'position',[tmp(1) tmp(2)-20 tmp(3) tmp(4)],...
                   'string','Left/Right',...
				   'backgroundcolor',GD.col.mainfig, ...
                   'value',0);

% File listbox
tmp = get(H.FLIPUIGRP,'position');
H.FILELBOX = uicontrol('parent',H.OPTUIPANEL,...
                       'units','pixel',...
                       'style','listbox',...
                       'position',[tmp(1)+tmp(3)+10 tmp(2) fig_w-tmp(1)-tmp(3)-35 ...
                   200],...
                       'backgroundcolor','w',...
                       'string',Dat.fullfiles,...
                       'Min',0,'Max',2);
tmp = get(H.FILELBOX,'position');
files_tx = uicontrol('parent',H.OPTUIPANEL,...
                     'units','pixel',...
                     'style','text',...
                     'position',[tmp(1) tmp(2)+tmp(4) 150 12],...
                     'string','Files (slices)',...
                     'horizontalalign','left',...
					 'backgroundcolor',GD.col.mainfig, ...
                     'fontweig','bold');

% Select all and select none buttons
tmp = get(H.FILELBOX,'position');
H.SELALL = uicontrol('parent',H.OPTUIPANEL,...
                     'units','pixel',...
                     'style','pushbutton',...
                     'position',[tmp(1) tmp(2)-25-5 tmp(3)/2-5 25],...
                     'string','Select All',...
                     'callback',@l_SelectAll);
tmp = get(H.SELALL,'position');
H.SELNONE = uicontrol('parent',H.OPTUIPANEL,...
                      'units','pixel',...
                      'style','pushbutton',...
                      'position',[tmp(1)+tmp(3)+10 tmp(2) tmp(3) 25],...
                      'string','Select None',...
                      'userdata',H.FILELBOX,...
                      'callback','set(get(gcbo,''userdata''),''value'',[])');

% OK button
H.OKBTN = uicontrol('parent',H.FH,...
                    'units','pixel',...
                    'position',[350 5 70 30],...
                    'string','OK',...
                    'callback',@l_OKCallBack);

% Cancel button
tmp = get(H.OKBTN,'position');
H.CANCELTN = uicontrol('parent',H.FH,...
                       'units','pixel',...
                       'position',[tmp(1)+tmp(3)+5 tmp(2) tmp(3) tmp(4)],...
                       'string','Cancel',...
                       'callback','delete(gcbf)');

end % function l_DrawGUI()


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Select All Files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_SelectAll(h,evd)
  lbox_length = length(get(H.FILELBOX,'string'));
  set(H.FILELBOX,'value',1:lbox_length);
end % function l_SelectAll(h,


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Selection Changed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_SelectionChanged(h,evd)

if h==H.ROTUIGRP
  
  if evd.NewValue==H.ROTNONE
    set([H.FLIPNONE,H.FLIPUD,H.FLIPLR],'enable','on')
  else
    set([H.FLIPNONE,H.FLIPUD,H.FLIPLR],'enable','off')
  end
  
elseif h==H.FLIPUIGRP
  
  if evd.NewValue==H.FLIPNONE
    if ~squareImages
      set([H.ROTNONE,H.ROT180],'enable','on')
    else
      set([H.ROTNONE,H.ROT90,H.ROT180,H.ROT270],'enable','on')
    end
  else
    set([H.ROTNONE,H.ROT90,H.ROT180,H.ROT270],'enable','off')
  end
  
end

end % function



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OK Selected
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l_OKCallBack(h,evd)

cancel = false;
val = get(H.FILELBOX,'value');
str = get(H.FILELBOX,'string');
RotSelObject = get(H.ROTUIGRP,'SelectedObject');
FlipSelObject = get(H.FLIPUIGRP,'SelectedObject');

%% Flipping vector:
%% 0 = none
%% 1 = Up/Down
%% 2 = Left/Right
Out.Flip = zeros(1,length(str));
if ~isempty(val)
  if FlipSelObject==H.FLIPUD
    Out.Flip(val)=1;
  elseif FlipSelObject==H.FLIPLR
    Out.Flip(val)=2;
  end
end

%% Rotate vector
%% 0 = none
%% 1 = 90
%% 2 = 180
%% 3 = 270
Out.Rotate = zeros(1,length(str));
if ~isempty(val)
  if RotSelObject==H.ROT90
    Out.Rotate(val)=1;
  elseif RotSelObject==H.ROT180
    Out.Rotate(val)=2;
  elseif RotSelObject==H.ROT270
    Out.Rotate(val)=3;
  end
end

%% Delete figure window
delete(H.FH)

end % function l_OKCallBack(h,

end
