function [coX,coiG,coFaces] = xDistributeIP(X,g,cg,overlap,bc,activeBnd)

%
% Distributes inner products stored in X as returned from
% xComputeMimeticIP so that basis functions can be constructed by calling 
% xGenerateCoarseSystem using coX, coiG and coFaces.
%
% PARAMETERS:
%
% - X -- Distributed cell array. Returned from xComputeMimeticIP. Contains
%     inverse ip_simple inner products for each cell.
%
% - G -- Composite variable. geometry structure from MRST.
%
% - CG -- composite variable. Coarse geometry structure from MRST.
%
% - Overlap -- Composite variable. Number of fine-grid cells in each
%     physical direction with which to extend the supporting domain 
%     of any given basis functions.
%
% - BC -- Composite variable. Boundary conditions structure from MRST.
%
% - activeBnd -- Composite variable. Vector of active coarse boundary
%      faces.
%      use activeBnd = [] (only no-flow BC) as default. We remark
%      that coarse faces with prescribed fine-scale BC's
%      are always considered active.
%
% RETURNS:
%
% coX -- Distributed cell array. Local part on each worker contains all
%   inner products that will be needed on that worker to construct basis
%   functions.
%
% coiG -- Distributed cell array. One cell pr. worker. that is, for worker 
%   w we have that coiG{w}{i} == k means that coX{w}{i} is inner product 
%   of global cell k.
%
% coFaces -- Distributed array. Face indices.
%   running: 
%     lpFaces = getLocalPart(coFaces);
%   on a worker, means that lpFaces(:) will be calculated on that worker.
%   (An that the needed ips are available when this function has been
%   run.)
%
% EXAMPLE:
%   See xExample for a complete example of usage. 
% 
% SEE ALSO:
%   xExample, xInitWorkers, xBroadcast, xComputeMimeticIP, 
%   xGenerateCoarseSystem, xEvalBasisFunc

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

spmd
  if labindex == 1    
    if any(activeBnd > cg.faces.num), % from MRST
      nonext = activeBnd(activeBnd > cg.faces.num);
      s = ''; if numel(nonext) > 1, s = 's'; end
      error(id('CoarseFace:NonExistent'), ...
              ['Cowardly refusing to assign basis function on ', ...
               'non-existent coarse face%s: [%s].'], s, int2str(nonext));  
      elseif any(~any(cg.faces.neighbors(activeBnd,:) == 0, 2)),
        error(id('ActiveBndFaces:NotBoundary'), ...
                ['I am confused. At least one of the purported ', ...
                 '''activeBndFaces'' isn''t. A boundary face, that is...']); 
    end
  end
  
  % used to form coFaces.
  if ~isempty(bc)
    activeBnd = unique([activeBnd; has_bc(g, cg, bc)]);
  end
  
  % distributor
  dist = codistributor1d(1,ones(numlabs,1)',[numlabs 1]);

  % note that activeBnd is a replicated (composite) array
  % otherwise distributing this way will not work.
  coFaces = codistributed([find(all(cg.faces.neighbors>0, 2));activeBnd]);
  lpFaces = getLocalPart(coFaces);  % global indices of faces on worker
  
  blks = unique(reshape(cg.faces.neighbors(lpFaces,:),[],1)); % ditto blks
  blks = blks(blks>0);
  lpiG = blk2cells(g,cg,overlap,blks); % ditto cells
  
  % distributed array of global cell indices
  coiG = codistributed.build({lpiG},dist,'noCommunication');

  % temporarily gather ips and cell indices to worker 1
  % note that gG{w} conains global indices of cells that will be needed on 
  % worker 'w'
  gX = gather(X,1); gG = gather(coiG,1);
  
  % distribute needed ips to each worker
  if labindex == 1,
    lpX = gX(gG{1});
    for lab = 2:numlabs
      labSend(gX(gG{lab}),lab);
    end
  else
    lpX = labReceive(1);
  end
  % store required ips in distributed array coX.
  % uses cell structure to avoid communicating the number of ips pr. worker.
  % (would be needed to construct a complete distributor.)
  % alternative is to let coX be a composite array.
  coX = codistributed.build({lpX}, dist,'noCommunication');
end
end

function res = blk2cells(g,cg,overlap,blk) % adapted from MRST (sub_cells)
% returns cell indices of blocks in blk

sub_c = sparse(1 : g.cells.num, cg.partition, 1,...
                   g.cells.num, max(cg.partition));  % == cg.cells.subCells
if overlap > 0,
  n = double(g.faces.neighbors(all(g.faces.neighbors > 0, 2), :));
  n = sparse([n(:,1); n(:,2); (1 : g.cells.num).'], ...
             [n(:,2); n(:,1); (1 : g.cells.num).'], 1, g.cells.num, g.cells.num);

  % BFS to discover immediate neighbours in overlap region.
  for o = 1 : overlap, sub_c = n * sub_c; end
end
sub_c = logical(sub_c);
res =  find(any(sub_c(:,blk), 2)); % get cells of block 'blk'
end

function ix = has_bc(g, cg, bc) % adapted from MRST
[nsub, sub] = subFaces(g, cg);
f2c         = sparse(sub, 1, rldecode((1 : double(cg.faces.num)) .', nsub));

% Note:
%   This code assumes that any Neumann conditions specified on the fine
%   scale grid amounts to flow only one way across a coarse face.  If the
%   Neumann conditions amount to zero net flow across a coarse face, the
%   resulting flux basis function is undefined...
%
flow_support                 = false([cg.faces.num, 1]);
flow_support(f2c([bc.face])) = true;

ix = find(flow_support);
end
