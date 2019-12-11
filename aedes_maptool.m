function aedes_maptool()
% AEDES_MAPTOOL - A GUI for calculating various parameter maps
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
Dat = [];

% Draw GUI
H = l_DrawGUI();


%%%%%%%%%%%%%%%%%%%%%%%
%% Draw GUI
%%%%%%%%%%%%%%%%%%%%%%%
function H=l_DrawGUI()

% Get screensize
scrsz = get(0,'screensize');
fig_w = 650;
fig_h = 580;
fig_pos = [scrsz(3)/2-fig_w/2 scrsz(4)/2-fig_h/2 ...
          fig_w fig_h];

%% Load default font and colors
% DefaultColor=get(0,'DefaultUicontrolBackgroundcolor');
GD=aedes_gui_defaults;
% GD.col.frame = DefaultColor;

% Main figure --------------------------
H.FIG = figure('units','pixels',...
  'position',fig_pos,...
  'NumberTitle','off',...
  'Name','Aedes Map Tool',...
  'Menubar','none',...
  'toolbar','none',...
  'color',GD.col.mainfig, ...
  'doublebuffer','on',...
  'Handlevisibility','off');

% UiToolbar ----------------------
try
  tmp=load('aedes_cdata.mat');
catch
  error('Could not open aedes_cdata.mat. Aborting...')
end
Dat.cdata = tmp.cdata;
hUiToolbar = uitoolbar('parent',H.FIG);
H.UITOOLBAR_OPEN = uipushtool('parent',hUiToolbar,...
  'cdata',Dat.cdata.open2,...
  'tooltip','Open options file',...
  'clickedcallback','');
H.UITOOLBAR_SAVE = uipushtool('parent',hUiToolbar,...
  'cdata',Dat.cdata.save2_small,...
  'tooltip','Save options',...
  'clickedcallback','');



%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%

% Options uipanel ----------------------
H.OPTIONS_UIPANEL = uipanel('parent',H.FIG,...
  'units','pixels',...
  'position',[5 40 fig_w-10 300],...
  'BackGroundColor',GD.col.frame,...
  'title','Options',...
  'fontweight','bold');
tmp = get(H.OPTIONS_UIPANEL,'position');

H.OUTDIR_UIPANEL = uipanel('parent',H.OPTIONS_UIPANEL,...
  'units','pixels',...
  'position',[5 5 tmp(3)-10 80],...
  'BackGroundColor',GD.col.frame);
tmp = get(H.OUTDIR_UIPANEL,'position');

% Output directory editbox ---------------

% Get default output folder from preferences
try
  DefaultOutputFolder = getpref('Aedes','MapToolDefaultOutputFolder');
catch
  if isunix
	DefaultOutputFolder = getenv('HOME');
  else
	DefaultOutputFolder = getenv('USERPROFILE');
  end
end
H.OUTDIR_EDIT = uicontrol('parent',H.OUTDIR_UIPANEL,...
  'style','edit',...
  'units','pixels',...
  'position',[5 5 tmp(3)-60 22],...
  'backgroundcolor',GD.col.edit,...
  'enable','off',...
  'horizontalalign','left',...
  'string',DefaultOutputFolder);
tmp = get(H.OUTDIR_EDIT,'position');

% Output directory browse btn ------------
H.OUTDIR_BROWSE = uicontrol('parent',H.OUTDIR_UIPANEL,...
  'style','pushbutton',...
  'units','pixels',...
  'position',[tmp(1)+tmp(3)+3 tmp(2) 47 tmp(4)],...
  'string','...',...
  'enable','off');
tmp = get(H.OUTDIR_EDIT,'position');

% Input dir as output dir
H.OUTDIR_RADIO1 = uicontrol('parent',H.OUTDIR_UIPANEL,...
  'units','pixels',...
  'position',[tmp(1) tmp(2)+tmp(4)+3 350 tmp(4)],...
  'style','radio',...
  'value',0,...
  'String','Use custom folder',...
  'BackGroundColor',GD.col.frame,...
  'callback',{@l_SelectOutputDir,'custom'});
tmp = get(H.OUTDIR_RADIO1,'position');

H.OUTDIR_RADIO2 = uicontrol('parent',H.OUTDIR_UIPANEL,...
  'units','pixels',...
  'position',[tmp(1) tmp(2)+tmp(4)+3 tmp(3) tmp(4)],...
  'style','radio',...
  'value',1,...
  'String','Use data input folder as output folder',...
  'BackGroundColor',GD.col.frame,...
  'callback',{@l_SelectOutputDir,'input'});

% Maptype and fit values uipanel ------------------
tmp = get(H.OUTDIR_UIPANEL,'position');
panel_w = 215;
panel_h = 190;
H.MAPTYPE_UIPANEL = uipanel('parent',H.OPTIONS_UIPANEL,...
  'units','pixels',...
  'position',[5 tmp(2)+tmp(4)+5 panel_w panel_h],...
  'backgroundcolor',GD.col.frame);
tmp = get(H.MAPTYPE_UIPANEL,'position');

% Maptype popup --------------------------
H.MAPTYPE_TX = uicontrol('parent',H.MAPTYPE_UIPANEL,...
  'units','pixels',...
  'position',[5 tmp(4)-10-15 70 15],...
  'style','text',...
  'String','Map type',...
  'backgroundcolor',GD.col.frame,...
  'horizontalalign','left');
tmp = get(H.MAPTYPE_TX,'position');
maptypes = l_GetMapTypes;
H.MAPTYPE_POPUP = uicontrol('parent',H.MAPTYPE_UIPANEL,...
  'units','pixels',...
  'position',[tmp(1) tmp(2)-20-3 200 20],...
  'style','popup',...
  'string',maptypes);
tmp = get(H.MAPTYPE_TX,'position');

% Fit values edit and popup -------------------
H.FITVALS_TX = uicontrol('parent',H.MAPTYPE_UIPANEL,...
  'units','pixels',...
  'position',[tmp(1) tmp(2)-15-40 tmp(3) 15],...
  'style','text',...
  'String','Fit values',...
  'backgroundcolor',GD.col.frame,...
  'horizontalalign','left');
tmp = get(H.FITVALS_TX,'position');

H.FITVALS_EDIT = uicontrol('parent',H.MAPTYPE_UIPANEL,...
  'units','pixels',...
  'position',[tmp(1) tmp(2)-22-3 200 22],...
  'style','edit',...
  'String','',...
  'backgroundcolor',GD.col.edit,...
  'horizontalalign','left');
tmp = get(H.FITVALS_EDIT,'position');

% % Fit value multiplication
% H.FITVALS_MULTI = uicontrol('parent',H.MAPTYPE_UIPANEL,...
%   'units','pixels',...
%   'position',[tmp(1) tmp(2)-22-2 tmp(3) 22],...
%   'style','popup',...
%   'String',{'Multiplication: 1',...
%   'Multiplication: 1000',...
%   'Multiplication: 100',...
%   'Multiplication: 10',...
%   'Multiplication: 0.1',...
%   'Multiplication: 0.01',...
%   'Multiplication: 0.001',...
%   'Multiplication: custom'});
% tmp = get(H.FITVALS_MULTI,'position');

% Use procpar field --------------------------
H.FITVALS_USE_PROCPAR = uicontrol('parent',H.MAPTYPE_UIPANEL,...
  'units','pixels',...
  'position',[tmp(1) tmp(2)-22-2 tmp(3) 22],...
  'style','popup',...
  'tooltip',sprintf('Select the procpar field\n that contains the fit values'),...
  'String',{'Use procpar field...'},...
  'horizontalalign','left');
tmp = get(H.FITVALS_USE_PROCPAR,'position');
H.FITVALS_USE_PROCPAR_SAVE = uicontrol('parent',H.MAPTYPE_UIPANEL,...
  'units','pixels',...
  'position',[tmp(1) tmp(2)-20-2 tmp(3) 20],...
  'style','checkbox',...
  'tooltip',sprintf('Use this procpar field as\n default for this map type'),...
  'String','Use as default for this map type',...
  'horizontalalign','left',...
  'backgroundcolor',GD.col.frame,...
  'fontsize',8);


%% Output file name and mask file controls
tmp=get(H.MAPTYPE_UIPANEL,'position');
H.OUTPUTFILE_UIPANEL = uipanel('parent',H.OPTIONS_UIPANEL,...
  'units','pixels',...
  'position',[tmp(1)+tmp(3)+5 tmp(2) 270 tmp(4)],...
  'backgroundcolor',GD.col.frame);
tmp=get(H.OUTPUTFILE_UIPANEL,'position');

% Output file name text and edit
H.OUTPUTFILE_TX = uicontrol('parent',H.OUTPUTFILE_UIPANEL,...
  'units','pixels',...
  'position',[5 tmp(4)-15-10 120 15],...
  'style','text',...
  'String','Output File Name',...
  'backgroundcolor',GD.col.frame,...
  'horizontalalign','left');
tmp = get(H.OUTPUTFILE_TX,'position');
outputfile_tooltip = ...
  sprintf(['The following special formatting can be\n',...
  'used in this dialog:\n',...
  '%%m = map type identifier (t2map, t1rmap, etc.)\n',...
  '%%f = file name from input folder (e.g. MyData.fid -> %%f = MyData)\n',...
  ' \n',...
  'Note that numbering and file extension are handled automatically.']);
H.OUTPUTFILE_EDIT = uicontrol('parent',H.OUTPUTFILE_UIPANEL,...
  'units','pixels',...
  'position',[tmp(1) tmp(2)-22-3 250 22],...
  'style','edit',...
  'String','%m_%f',...
  'backgroundcolor',GD.col.edit,...
  'horizontalalign','left',...
  'tooltip',outputfile_tooltip);
tmp = get(H.OUTPUTFILE_EDIT,'position');

% Overwrite without warning checkbox
H.OUTPUTFILE_OVERWRITE = uicontrol('parent',H.OUTPUTFILE_UIPANEL,...
  'units','pixels',...
  'position',[tmp(1) tmp(2)-20-2 tmp(3) 20],...
  'style','checkbox',...
  'value',0,...
  'String','Overwrite without warning',...
  'backgroundcolor',GD.col.frame,...
  'horizontalalign','left');
tmp = get(H.OUTPUTFILE_OVERWRITE,'position');

% Fit type text
H.FITTYPE_TX = uicontrol('parent',H.OUTPUTFILE_UIPANEL,...
  'units','pixels',...
  'position',[5 tmp(2)-15-20 120 15],...
  'style','text',...
  'String','Fit type',...
  'backgroundcolor',GD.col.frame,...
  'horizontalalign','left');
tmp = get(H.FITTYPE_TX,'position');

% Linearized or non-linear fit
H.LINEARFIT_RADIO = uicontrol('parent',H.OUTPUTFILE_UIPANEL,...
  'units','pixels',...
  'position',[tmp(1) tmp(2)-22-3 200 22],...
  'style','radio',...
  'value',1,...
  'string','Linearized fit',...
  'backgroundcolor',GD.col.frame);
tmp = get(H.LINEARFIT_RADIO,'position');
H.NONLINEARFIT_RADIO = uicontrol('parent',H.OUTPUTFILE_UIPANEL,...
  'units','pixels',...
  'position',[tmp(1) tmp(2)-22-3 tmp(3) 22],...
  'style','radio',...
  'value',0,...
  'string','Non-linear fit',...
  'backgroundcolor',GD.col.frame);
tmp = get(H.NONLINEARFIT_RADIO,'position');

% % Initial values -------------------------------
% H.INITIALVAL_TX = uicontrol('parent',H.OUTPUTFILE_UIPANEL,...
%   'units','pixels',...
%   'position',[35 tmp(2)-15-3 90 15],...
%   'style','text',...
%   'String','Initial values',...
%   'backgroundcolor',GD.col.frame,...
%   'horizontalalign','left');
% tmp = get(H.INITIALVAL_TX,'position');
% 
% H.INITIALVAL_POPUP = uicontrol('parent',H.OUTPUTFILE_UIPANEL,...
%   'units','pixels',...
%   'position',[tmp(1)+tmp(3) tmp(2)-2 70 22],...
%   'style','popup',...
%   'string',{'Auto','Custom'},...
%   'enable','off');
% tmp = get(H.INITIALVAL_POPUP,'position');
% H.INITIALVAL_EDIT = uicontrol('parent',H.OUTPUTFILE_UIPANEL,...
%   'units','pixels',...
%   'position',[tmp(1)+tmp(3) tmp(2) 60 22],...
%   'style','edit',...
%   'String','',...
%   'backgroundcolor',GD.col.edit,...
%   'horizontalalign','left',...
%   'enable','off');

% % Mask File ---------------------------------
% H.MASKFILE_TX = uicontrol('parent',H.MASK_UIPANEL,...
%   'units','pixels',...
%   'position',[tmp(1) tmp(2)-15-20 tmp(3) 15],...
%   'style','text',...
%   'String','Mask File',...
%   'backgroundcolor','r',...GD.col.frame,...
%   'horizontalalign','left');
% tmp = get(H.MASKFILE_TX,'position');
% H.MASKFILE_EDIT = uicontrol('parent',H.MASK_UIPANEL,...
%   'units','pixels',...
%   'position',[tmp(1) tmp(2)-22 tmp(3) 22],...
%   'style','edit',...
%   'String','',...
%   'backgroundcolor',GD.col.edit,...
%   'horizontalalign','left');

% Add options to selected files -------------
tmp = get(H.OUTPUTFILE_UIPANEL,'position');
btn_h = 30;
btn_w = 130;
H.ADD_OPT_TO_SEL = uicontrol('parent',H.OPTIONS_UIPANEL,...
  'units','pixels',...
  'position',[tmp(1)+tmp(3)+5 tmp(2)+tmp(4)-btn_h btn_w btn_h],...
  'style','pushbutton',...
  'string','Add to selected',...
  'tooltip','Add options to selected files');
tmp = get(H.ADD_OPT_TO_SEL,'position');
H.ADD_OPT_TO_ALL = uicontrol('parent',H.OPTIONS_UIPANEL,...
  'units','pixels',...
  'position',[tmp(1) tmp(2)-btn_h-3 tmp(3) tmp(4)],...
  'style','pushbutton',...
  'string','Add to all',...
  'tooltip','Add options to all files');


%% FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Files uipanel ---------------------
tmp = get(H.OPTIONS_UIPANEL,'position');
H.FILES_UIPANEL = uipanel('parent',H.FIG,...
  'units','pixels',...
  'position',[tmp(1) tmp(2)+tmp(4)+5 fig_w-10 220],...
  'BackGroundColor',GD.col.frame,...
  'title','Data Files',...
  'fontweight','bold');
tmp = get(H.FILES_UIPANEL,'position');

% Quit button --------------------------
btn_h = 30;
btn_w = 90;
H.QUIT_BTN = uicontrol('parent',H.FIG,...
  'units','pixels',...
  'position',[fig_w-btn_w-5 5 btn_w btn_h],...
  'style','pushbutton',...
  'string','Quit');
tmp = get(H.QUIT_BTN,'position');

% Fit maps btn-------------------------
H.FIT_MAPS_BTN = uicontrol('parent',H.FIG,...
  'units','pixels',...
  'position',[tmp(1)-tmp(3)-3 tmp(2) tmp(3) tmp(4)],...
  'style','pushbutton',...
  'string','Fit Maps');

% Files listbox -------------------------
tmp = get(H.FILES_UIPANEL,'position');
H.FILES_LBOX = uicontrol('parent',H.FILES_UIPANEL,...
  'style','listbox',...
  'string',{'*/home/tjniskan/MyMapData/070907/MapFromARat.fid/fid',...
  '*/home/tjniskan/MyMapData/070907/MapFromARat.fid/fid',...
  '*/home/tjniskan/MyMapData/090907/MapFromARat.fid/fid',...
  '*/home/tjniskan/MyMapData/100907/MapFromARat.fid/fid',...
  '*/home/tjniskan/MyMapData/120907/MapFromARat.fid/fid',...
  '/home/tjniskan/MyMapData/150907/MapFromARat.fid/fid',...
  '/home/tjniskan/MyMapData/180907/MapFromARat.fid/fid',...
  '/home/tjniskan/MyMapData/230907/MapFromARat.fid/fid',...
  '/home/tjniskan/MyMapData/011007/MapFromARat.fid/fid',...
  '/home/tjniskan/MyMapData/051007/MapFromARat.fid/fid',...
  '/home/tjniskan/MyMapData/131007/MapFromARat.fid/fid'},...
  'units','pixel',...
  'position',[5 5 tmp(3)-10 tmp(4)-70],...
  'backgroundcolor',GD.col.listbox,...
  'Max',3,'Min',1);
tmp = get(H.FILES_LBOX,'position');

% Add files btn --------------------------
H.FILES_ADD = uicontrol('parent',H.FILES_UIPANEL,...
  'style','pushbutton',...
  'units','pixels',...
  'position',[5 tmp(2)+tmp(4)+3 40 40],...
  'tooltip','Add files to list',...
  'cdata',Dat.cdata.add);
tmp = get(H.FILES_ADD,'position');

% Remove selected ------------------------
H.FILES_REMSEL = uicontrol('parent',H.FILES_UIPANEL,...
  'style','pushbutton',...
  'units','pixels',...
  'position',[tmp(1)+tmp(3)+3 tmp(2) tmp(3) tmp(4)],...
  'tooltip','Remove selected files from list',...
  'cdata',Dat.cdata.delete);
tmp = get(H.FILES_REMSEL,'position');

% Remove all -----------------------------
H.FILES_REMALL = uicontrol('parent',H.FILES_UIPANEL,...
  'style','pushbutton',...
  'units','pixels',...
  'position',[tmp(1)+tmp(3)+3 tmp(2) tmp(3) tmp(4)],...
  'tooltip','Remove all files from list',...
  'cdata',Dat.cdata.deleteall);
tmp = get(H.FILES_REMALL,'position');




end % function l_DrawGUI()

%%%%%%%%%%%%%%%%%%%%%
%% Get Map Types
%%%%%%%%%%%%%%%%%%%%%
  function maptypes = l_GetMapTypes()
	maptypes = {'T2','T1','T2 rho','T1 rho','ADC',...
	  'Perfusion','B1'};
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Callback for output folder selection radiobuttons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function l_SelectOutputDir(h,evd,opt)
	
	if strcmpi(opt,'custom')
	  set(H.OUTDIR_EDIT,'enable','on')
	  set(H.OUTDIR_BROWSE,'enable','on')
	  set(H.OUTDIR_RADIO1,'value',1)
	  set(H.OUTDIR_RADIO2,'value',0)
	elseif strcmpi(opt,'input')
	  set(H.OUTDIR_EDIT,'enable','off')
	  set(H.OUTDIR_BROWSE,'enable','off')
	  set(H.OUTDIR_RADIO1,'value',0)
	  set(H.OUTDIR_RADIO2,'value',1)
	end
	
  end

end % aedes_maptool()
