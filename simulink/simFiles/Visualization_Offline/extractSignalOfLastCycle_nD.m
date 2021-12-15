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
% :Adapted from: Dylan Eijkelhof (d.eijkelhof@tudelft.nl)

function y = extractSignalOfLastCycle_nD( signal, sample_count_last_cycle,  simInit, number_of_clycles)
%Extract n-D data from last converged power cycle.
%
% :param signal: Full n-D timeseries simulation output,size= nx1xtimeseries
% :param sample_count_last_cycle: Simulation timeseries 'sample_count_last_cycle'
% :param simInit: Simulation input structure 'simInit'
%
% :returns:
%           - **y** - Extracted power cycle n-D timeseries.
%
%------------- BEGIN CODE --------------
if nargin < 4
    number_of_clycles = 1;
end
time_window_last_cylce = sum(sample_count_last_cycle(end:-1:end-number_of_clycles+1)) * simInit.Ts_power_conv_check;
idx_time_window_start = find(signal.Time >= signal.Time(end)-time_window_last_cylce,1);

s = size(signal.Data);
if numel(s)==3
    one_pos = find(s==1,1);
    snew = s;
    snew(one_pos)=[];
    i_pos = find(s==min(snew));
    i_pos = i_pos(end);
    for i=1:s(i_pos)
        if i_pos==1
            if one_pos == 2
                y.Data(i,:) = signal.Data(i,1,idx_time_window_start:end);
            else
                y.Data(i,:) = signal.Data(i,idx_time_window_start:end,1);
            end
        elseif i_pos==2
            if one_pos == 1
                y.Data(i,:) = signal.Data(1,i,idx_time_window_start:end);
            else
                y.Data(i,:) = signal.Data(idx_time_window_start:end,i,1);
            end
        else
            if one_pos == 1
                y.Data(i,:) = signal.Data(1,idx_time_window_start:end,i);
            else
                y.Data(i,:) = signal.Data(idx_time_window_start:end,1,i);
            end
        end
    end
%     for i=1:size(signal.Data,1)
%         y.Data(i,:) = signal.Data(i,1,idx_time_window_start:end);
%     end
elseif numel(s)==2
    [~,i_pos] = min(s);
    for i=1:s(i_pos)
        if i_pos==1
            y.Data(i,:) = signal.Data(i,idx_time_window_start:end);
        else
            y.Data(i,:) = signal.Data(idx_time_window_start:end,i)';
        end
    end
else
    error('Input is not compatible')
end
    
y.Time = signal.Time(idx_time_window_start:end);
%------------- END CODE --------------
end