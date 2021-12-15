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
% :Adapted from: Sebastian Rapp, TU Delft

function makevideo(output_name, preamble)

    if nargin <2
        preamble = 'video/vidPic_';
    end
    if nargin < 1
        output_name = 'pumping_cycle';
    end
       
    writerObj = VideoWriter(output_name,'MPEG-4');

    open(writerObj);
    
    K = 1;
    datenum_prev = 0;
    
    while K <= inf
        filename = sprintf([preamble,'%d.png'], K);
        if ~isfile(filename)
            break
        end
        FileInfo = dir(filename);
        if FileInfo.datenum < datenum_prev 
            break
        else
            datenum_prev = FileInfo.datenum;
        end
        thisimage = imread(filename);
        writeVideo(writerObj, thisimage);
        K = K +1; 
    end
    close(writerObj);
end