function coX = xComputeMimeticIP(g,rock,facetrans)
% 
% calculates inverse ip_simple inner products in parallel and stores in 
% distributed cell array coX.
% that is, the ip of cell i is stored in coX{i}. 
% 
% PARAMETERS:
%
% - g -- composite variable. geometry structure as in MRST.
%
% - rock -- composite variable. rock structure as in MRST. 
%
% - facetrans -- composite variable. 
%     If facetrans is specified, the innerproduct is modified to account
%     for a face transmissibilities specified as
%       [faces, face_transmissibilities]
%   use facetrans = [], if there are no face transmissibilities.
% 
% RETURNS:
%
% - coX -- distributed cell array. coX{i} contains inverse ip_simple inner
%     product of cell i.    
%     moreover, to get the full BI from coX, run the following code on
%     host:
%       dimProd = double(diff(g.cells.facePos));
%       X = gather(coX);
%       [ind1, ind2] = blockDiagIndex(dimProd, dimProd);
%       n = size(G.cells.faces, 1);
%       BI = sparse(ind1,ind2,vertcat(res_invip{:}),n,n);
%
%
% use xDistributeIP to distribute inner products to active workers in order
% to construct basis functions in parallel with xGenerateCoarseSystem.
%

spmd
  % used on all workers
  dim = size(g.nodes.coords, 2);
  % cell faces
  cf = double(g.cells.faces(:,1));
  % cell number of face
  cellno = rldecode(1 : g.cells.num, diff(g.cells.facePos), 2) .';
  % direction of face vector (in/out)
  sgn = 2*double( g.faces.neighbors(cf, 1) == cellno) - 1;

  % run in parallel on all workers. reduces communication.
  areas      = g.faces.areas(cf);
  centVecs   = g.faces.centroids(cf,:)-g.cells.centroids(cellno,:);
  outNormals = bsxfun(@times,g.faces.normals(cf,:),sgn);
  perm       = permTensor(rock, dim);

  assert(size(perm,1) == g.cells.num, ...
         ['Permeability must be defined in active cells only.\n', ...
          'Got %d tensors, expected %d (== number of cells).'],   ...
          size(perm,1), g.cells.num);

  % Compute half-face transmissibilities from face transmissibilities.
  if ~isempty(facetrans),
     if ~all(facetrans(:,2) > 0),
        disp('\n\tAll face transmissibilities should be positive.\n\t\t\t');
     end
     hft = zeros(g.faces.num, 1); % distributed to all workers later
     hft(facetrans(:,1)) = facetrans(:,2);
     hft = hft.*sum(g.faces.neighbors ~= 0, 2);
     i=hft~=0; hft(i) = 1./hft(i);
     hft = hft(g.cells.faces(:,1));
  else
     hft = zeros(size(g.cells.faces(:,1)));
  end

  dimProd = double(diff(g.cells.facePos));
  cumProd = [0;cumsum(dimProd)];

  % simplest (?) way to naively distribute cell indices 
  % codistributed(a); a must be replicated array
  coiG = codistributed( (1:g.cells.num)' );
  dist = getCodistributor(coiG); % distributor for building result
  lpiG = getLocalPart(coiG); % cell indices of cells on worker

  niG = numel(lpiG); % number of cells on worker
  lpX = cell(niG,1); % for storing ip on each worker
  
  % calculate inverse ip_simple for each cell on worker
  for i = 1:niG
    k = lpiG(i);
    nF = dimProd(k); ix = cumProd(k);
    pR = ix + 1 : ix + nF; % half-face indices

    K = reshape(perm(k,:), [dim, dim]);
    v = g.cells.volumes(k);

    a = areas(pR);
    C = centVecs(pR,:);
    N = outNormals(pR,:);
    
    this_ip = inv_ip_simple(a, v, K, C, N); % get inverse ip_simple of cell
    if any(hft(pR)),
      lpX{i,1} = reshape(inv(inv(this_ip) + diag(hft(pR))),[],1);
    else
      lpX{i,1} = reshape(this_ip,[],1);
    end
  end
  % store ip in distributed array lpX
  coX  = codistributed.build(lpX, dist,'noCommunication');
end
end

function W = inv_ip_simple(a, v, K, C, N) % adapted from MRST
t = 6 * sum(diag(K)) / size(K,2);
Q = orth(bsxfun(@times, C, a));
U = eye(length(a)) - Q*Q';
d = diag(a);
W = (N*K*N' + t*(d * U * d)) ./ v;
end
