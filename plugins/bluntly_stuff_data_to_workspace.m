function blunt_stuffer(DATA,ROI,AddInfo)
% Bluntly stuff the data to workspace, nevermind about existing variables
% with same name...
%
% This function is written for Aedes - A graphical tool for analyzing 
% medical images
%
% Copyright (C) 2011 mikko nissi <nissi@cmrr.umn.edu>
% 
% This program may be used under the terms of the GNU General Public
% License version 2.0 as published by the Free Software Foundation
% and appearing in the file LICENSE.TXT included in the packaging of
% this program.
%
% This program is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
% WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.


% Assign variables into workspace
assignin('base','DATA',DATA);
assignin('base','ROI',ROI);
assignin('base','AddInfo',AddInfo);


end
