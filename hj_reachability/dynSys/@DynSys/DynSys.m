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
% :Author: Mo Chen (mochen@sfu.ca)

classdef DynSys < handle
  % Dynamical Systems class; inherits from handle to make objects behave like
  % pointers
    
  properties % For bookkeepping and plotting
    nx          % Number of state dimensions
    nu          % Number of control inputs
    nd          % Number of disturbance dimensions
    
    x           % State
    u           % Recent control signal
    
    xhist       % History of state
    uhist       % History of control
    
    pdim        % position dimensions
    vdim        % velocity dimensions
    hdim        % heading dimensions
    
    %% Figure handles
    hpxpy           % Position
    hpxpyhist       % Position history
    hvxvy           % Velocity
    hvxvyhist       % Velocity history
    
    % Position velocity (so far only used in DoubleInt)
    hpv = cell(2,1);
    hpvhist = cell(2,1);
    
    % Data (any data that one may want to store for convenience)
    data
  end % end properties

  % No constructor in DynSys class. Use constructors in the subclasses
  
end % end class