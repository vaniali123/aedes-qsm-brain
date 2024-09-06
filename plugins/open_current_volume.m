function open_current_volume(DATA,ROI,AddInfo)
% This Aedes plugin calculates T1 map

% for some reason requires roi
% This function is written for Aedes
%
% copyright (c) 2008 Mikko Nissi <mikko.nissi@iki.fi>
%
%
% This program may be used under the terms of the GNU General Public
% License version 2.0 as published by the Free Software Foundation
% and appearing in the file LICENSE.TXT included in the packaging of
% this program.
%
% This program is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
% WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
% keyboard
try 
    % ask roilabels:
%     labels=aedes_inputdlg('surf1 surf2 bottom1 bottom2','Give roi labels:','s1 s2 b1 b2');
%     widh=aedes_inputdlg('Width of profile in pixels?','Give roi width:','7');
%     dim=aedes_inputdlg('dimension,x=1,y=2,z=3','Give correct number:','2');
%     corns=regexp(labels,'\w*','match');
%       keyboard
%     width=str2num(widh{1});
%     dimension = str2num(dim{1});
    %    keyboard;

    % Zero Infs and NaNs
    DATA{1}.FTDATA(isnan(DATA{1}.FTDATA))=0;
    DATA{1}.FTDATA(isinf(DATA{1}.FTDATA))=0;
    

   DATA{1}.FTDATA = DATA{1}.FTDATA(:,:,:,AddInfo.CurrentVol);

    
    aedes(DATA)
    
    

%    keyboard;
catch
  errordlg({'Could not do something. The following error was returned',...
           '',lasterr},'modal')
end