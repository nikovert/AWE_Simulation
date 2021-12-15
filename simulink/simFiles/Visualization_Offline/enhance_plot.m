% Copyright (C) 2021  Nikolaus Vertovec
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
% :Revision: 14-December-2021
% :Author: Nikolaus Vertovec (nikolaus.vertovec@eng.ox.ac.uk)
% :Adapted from: Dylan Eijkelhof (d.eijkelhof@tudelft.nl) and J. Nelson

function enhance_plot(fontname,fontsize,linewid,markersiz,lgd)
%Function to enhance MATLAB's lousy text choices on plots.  Sets the
%   current figure's Xlabel, Ylabel, Title, and all Text on plots, plus
%   the axes-labels to the "fontname" and "fontsize" input here where
%   the defaults have been set to 'times' and 16.
%   Also sets all plotted lines to "linewid" and all markers to size
%   "markersiz".  The defaults are 2 and 8.
%
% :param fontname:   (Optional,DEF='TIMES') FontName string to use MATLAB's ugly default is 'Helvetica'
% :param fontsize:   (Optional,DEF=16) FontSize integer to use MATLAB's tiny default is 10
% :param linewid:    (Optional,DEF=2) LineWidth integer to use MATLAB's skinny default is 0.5
% :param markersiz:  (Optional,DEF=8) MarkerSize integer to use MATLAB's squinty default is 6
% :param lgd:        
%             - (Optional, DEF=0)
%             - if is 0, doesn't change the legend
%             - if is 1, changes only the lines on the legend
%             - if is 2, changes both the lines and the text
%             - if is 3, changes only the text for all inputs, 
%             - if pass 0, use default
%             - if pass -1, use MATLAB's default
% :returns: 
%           - **vec_Abar** - Vector in rotated aerodynamic reference frame.

if (~exist('fontname','var')||(all(fontname==0) && isnumeric(fontname)))
  fontname = 'times';
elseif (fontname==-1)
  fontname = 'helvetica';
end
if (~exist('fontsize','var')||(fontsize==0))
  fontsize = 16;
elseif (fontsize==-1)
  fontsize=10;
end
if (~exist('linewid','var')||(linewid==0))
  linewid=2;
elseif (linewid==-1)
  linewid=0.5;
end
if (~exist('markersiz','var')||(markersiz==0))
  markersiz = 8;
elseif (markersiz==-1)
  markersiz = 6;
end
if (~exist('lgd','var')||(lgd==0)||(lgd<=-1))
  lgd=0;
end

box on; grid on;
Hf=gcf;
Ha=gca;
Hx=get(Ha,'XLabel');
Hy=get(Ha,'YLabel');
Ht=get(Ha,'Title');
set(Ha,'LineWidth',.75);
set(Hx,'fontname',fontname);
set(Hx,'fontsize',fontsize);
set(Hy,'fontname',fontname);
set(Hy,'fontsize',fontsize);
set(Ha,'fontname',fontname);
set(Ha,'fontsize',fontsize);
%set(Ha,'YaxisLocation','right')
%set(Ha,'YaxisLocation','left')
set(Ht,'fontname',fontname);
set(Ht,'fontsize',fontsize);
set(Hy,'VerticalAlignment','bottom');
set(Hx,'VerticalAlignment','cap');
set(Ht,'VerticalAlignment','baseline');
Hn = get(Ha,'Children');
n = length(Hn);
if n > 0
  typ = get(Hn,'Type');
  for j = 1:n
    if strcmp('text',typ(j,:))
      set(Hn(j),'fontname',fontname);
      set(Hn(j),'fontsize',fontsize);
    end
    if strcmp('line',typ(j,:))
      set(Hn(j),'LineWidth',linewid);
      set(Hn(j),'MarkerSize',markersiz);
    end
  end
end
%           legend:     (Optional, DEF=0) if is 0, doesn't change the legend
%                       if is 1, changes only the lines on the legend
%                       if is 2, changes both the lines and the text
%                       if is 3, changes only the text
if (lgd~=0)
  legh=legend;
  Hn=get(legh,'Children');
  n = length(Hn);
  if n > 0
    typ = get(Hn,'Type');
    for j = 1:n
      if (strcmp('text',typ(j,:)) && ((lgd==2)||(lgd==3)))
        set(Hn(j),'fontname',fontname);
        set(Hn(j),'fontsize',fontsize-2);
      end
      if (strcmp('line',typ(j,:)) && ((lgd==1)||(lgd==2)))
        set(Hn(j),'LineWidth',linewid);
        set(Hn(j),'MarkerSize',markersiz);
      end
    end
  end
end

  
  
  
figure(Hf);