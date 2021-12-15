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
% :Adapted from: Dylan Eijkelhof (d.eijkelhof@tudelft.nl) and Sebastian Rapp (s.rapp@tudelft.nl)

function y = extractSignalOfLastCycle2(signal, sample_count_last_cycle,  simInit, number_of_clycles)
%Extract 1D data from last converged power cycle.
%
% :param signal: Full 1D timeseries simulation output.
% :param sample_count_last_cycle: Simulation timeseries 'sample_count_last_cycle'.
% :param simInit: Simulation input structure 'simInit'.
%
% :returns:
%           - **y** - Extracted power cycle 1D timeseries.

%------------- BEGIN CODE --------------
if nargin < 4
    number_of_clycles = 1;
end

time_window_last_cylce = sum(sample_count_last_cycle(end:-1:end-number_of_clycles+1)) * simInit.Ts_power_conv_check;
idx_time_window_start = find(signal.Time >= signal.Time(end)-time_window_last_cylce,1);
y.Data = signal.Data(idx_time_window_start:end,:);
y.Time = signal.Time(idx_time_window_start:end);

%------------- END CODE --------------
end