function [g,rock,cg,mob,bc,src,overlap,facetrans,...
          weighting,activeBnd] = xInitWorkers(p)

% Initializes PCT workers and composite variables for use with the 
% xmsmfem module.
% All returned variables are Composite
% The following returned variables are assigned default values:
%   overlap   = 0;
%   bc        = [];
%   src       = [];
%   facetrans = zeros(0,2);
%   activeBnd = [];
%   weighting = 'perm';
% the remaining variables are empty.
% they need to be initialized. afterwards xBroadcast needs to be called.
% see xExample
%
% PARAMETERS:
% - p -- Number of wanted workers. Must be larger than 0. 
%     If an incorrect number of workers are already active, these are closed 
%     first. Fails if p workers can not be started.
%
% RETURNS:
% - G -- Composite variable. Grid structure from MRST.
%
% - rock -- Composite variable. Rock structure from MRST.
% 
% - CG -- composite variable. Coarse grid structure from MRST.
%
% - mob -- Composite variable. Total mobility. One scalar value for each
%     cell in the underlying (fine) model
%
% - bc -- Composite variable. Boundary condition structure from MRST.
%
% - src -- Composite variable. Explicit source contributions as defined by
%     'addSource'. 
%
% - overlap -- Composite variable. Number of fine-grid cells in each 
%     physical direction with which to extend the supporting domain of any
%     given basis functions.
%
% - facetrans -- composite variable. 
%     If facetrans is specified, the innerproduct is modified to account
%     for a face transmissibilities specified as
%       [faces, face_transmissibilities]
%
% - weighting -- Composite variable. Basis function driving source term as
%     supported by function 'evalBasisSource'.
%     use weighting = 'perm' as default value.
%
% - activeBnd -- Composite variable. Vector of active coarse boundary
%      faces.
%
% EXAMPLE:
%   See xExample for a complete example of usage. 
% 
% SEE ALSO:
%   xExample, xBroadcast, xComputeMimeticIP, 
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


assert(p > 0,'need to have at least one active worker.');

if matlabpool('size')~=p,
  if matlabpool('size'),matlabpool('close');end
  matlabpool('open',p);
end

assert(matlabpool('size')==p,'incorrect poolsize');

g         = Composite(); g(:)         = cell(1,p);
rock      = Composite(); rock(:)      = cell(1,p);
cg        = Composite(); cg(:)        = cell(1,p);
mob       = Composite(); mob(:)       = cell(1,p);
bc        = Composite(); bc(:)        = cell(1,p);
src       = Composite(); src(:)       = cell(1,p);
overlap   = Composite(); overlap(:)   = cell(1,p);
facetrans = Composite(); facetrans(:) = cell(1,p);
activeBnd = Composite(); activeBnd(:) = cell(1,p);
weighting = Composite(); weighting(:) = cell(1,p);

spmd
  if labindex == 1
    % initialize default values for Composite variables.
    % remember to call xBroadcast when you have finished setting geometry
    % and rock properties.
    % you will also need to call xBroadcast if you change any of the
    % following variables as well.
    % any changes need to be done on worker 1.
    overlap = 0;
    bc = [];
    src = [];
    facetrans = zeros(0,2);
    activeBnd = [];
    weighting = 'perm' ;
  end
end
end
