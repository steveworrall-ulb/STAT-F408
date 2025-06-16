
function [y transpose] = column(x)

% ThreshLab/General/column: makes sure that a vector has 1 column
%  Usage
%    y = column(x);
%    [y transpose] = column(x);
%  Inputs
%    x : row or column vector
%  Outputs
%    y : column vector
%    transpose : boolean/logical/binary variable transpose=(y==x')
%  Description
%  Note
%  See also
%    help row

if (size(x,1) < size(x,2) )
 y = x'; transpose = true;
else
 y = x; transpose = false;
end

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
