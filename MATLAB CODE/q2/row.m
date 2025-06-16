
function [y transpose] = row(x);

% ThreshLab/General/row: makes sure that a vector has 1 row
%  Usage
%    y = row(x);
%    [y transpose] = row(x);
%  Inputs
%    x : row or column vector
%  Outputs
%    y : row vector
%    transpose : boolean/logical/binary variable transpose=(y==x')
%  Description
%  Note
%  See also
%    help column


if (size(x,1) > size(x,2) )
 y = x'; transpose = true;
else
 y = x; transpose = false;
end

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
