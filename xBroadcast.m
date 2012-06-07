function [g,rock,cg,mob,bc,...
         src,overlap,facetrans,...
         weighting,activeBnd] = ...
         xBroadcast(g,rock,cg,mob,bc,...
         src,overlap,facetrans,...
         weighting,activeBnd)

% Helper function for the xmsmfem module for MRST.
% distributes initialized Composite variables before they can be used on workers
% See xExample to see how this function is used.
% 
% EXAMPLE:
%   See xExample for a complete example of usage. 
% 
% SEE ALSO:
%   xExample, xInitWorkers, xComputeMimeticIP, 
%   xDistributeIP, xGenerateCoarseSystem, xEvalBasisFunc

%{
A part of the xmsmfem module for MRST:
http://www.sintef.no/Projectweb/MRST/
Adapted from the msmfem module with the Parallel Computing Toolbox

Released under the GNU General Public License:
http://www.gnu.org/licenses/gpl.html
 
Written by
Anders Hoff 2012
http://master.andershoff.net
%}


g         = labBroadcast(1,g);
rock      = labBroadcast(1,rock);
cg        = labBroadcast(1,cg);
mob       = labBroadcast(1,mob);
bc        = labBroadcast(1,bc);
src       = labBroadcast(1,src);
overlap   = labBroadcast(1,overlap);
facetrans = labBroadcast(1,facetrans);
weighting = labBroadcast(1,weighting);
activeBnd = labBroadcast(1,activeBnd);

end
