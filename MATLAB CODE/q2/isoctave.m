
function B = isoctave();

% ThreshLab/General/isoctave: TRUE if called when running octave
%  Usage
%    B = isoctave;
%  Inputs
%    (none)
%  Outputs
%    B binary value (TRUE or FALSE, 1 or 0)
%  Description
%  Note
%  See also

B = exist('OCTAVE_VERSION','builtin') ~= 0;
