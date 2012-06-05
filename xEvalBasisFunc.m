function [coV, coP] = xEvalBasisFunc(g,cg,coX,coiG,...
                       coFaces,weight,mob,src,bc,overlap)
%
%
%
%
%
%

spmd
  linSolve = @mldivide;
  nBlk = max(cg.partition);
  avgmob = average_mobility(cg.partition,mob,g.cells.volumes);
  
  theta = get_coarse_weighting(g, cg.partition,weight,src);
  
  [f_bc, h_bc, is_dirichlet] = expand_bc(g,bc); 
  [sub_map, ncells, nhf_blk] = mappings(g,cg,mob,bc,overlap);

  dist     = getCodistributor(coFaces); % distributor
  lpFaces  = getLocalPart(coFaces); % indexes of local faces
  lnF      = length(lpFaces); % number of local faces on the worker
  press_ok = false(lnF,1); % to test for orthogonality
  lpV      = cell(lnF ,1); % store local basis of V
  lpP      = cell(lnF ,1); % store local basis of P
  numfaces = diff(g.cells.facePos);
  lpX      = getLocalPart(coX); % local ips
  lpiG     = getLocalPart(coiG); % local cell indices
  
  for f = 1:lnF
    face = lpFaces(f);
    blk  = cg.faces.neighbors(face,:); % Blocks connected to 'face'.
    blk  = blk(blk > 0);
    nblk = numel(blk);       % Number of blocks sharing 'face'.
    iG = sub_map.cells(blk); % Fine-scale cells      present in 'blk'
    iF = sub_map.hf(iG) ;    % Fine-scale half-faces present in 'blk'
    iH = sub_map.faces(iF) ; % Fine-scale faces      present in 'blk'

    niF = numel(iF); niG = numel(iG);

    % Extract linsys matrix components for sub-(hf,cells,faces).
    % as oppsed to evalBasisFunc theses matrices are constructed
    % locally to avid storing and distributing them.
    
    % note that lpX is a distributed cell array of inner products.
    [ind1, ind2] = blockDiagIndex(numfaces(iG));
    sBI = sparse(ind1, ind2, vertcat(lpX{1}{ismembc(lpiG{1},iG)}));

    n  = diff([g.cells.facePos(iG), g.cells.facePos(iG+1)], [], 2);
    sC = sparse(1:numel(iF), rldecode(1:numel(n), n, 2), 1);

    renum_f     = zeros([g.faces.num, 1]);
    renum_f(iH) = 1:numel(iH);
    sD = sparse(1:numel(iF), renum_f(g.cells.faces(iF,1)), 1);

    % Build linsys right hand side components for hybrid system.

    src_mult = zeros([nBlk,1]); sgn = [1,-1];
    src_mult(blk) = sgn(1:nblk);
    sG            = theta(iG) .* src_mult(cg.partition(iG));

    % 2) Define trivial vectors [f,h] (no external forces).
    %    Will be non-trivial in case of boundary conditions (below).
    sF  = zeros([niF, 1]); sH = zeros([numel(iH), 1]);

    % Update linsys components for presence of several phases...
    dmob = sub_map.dmob(iF, niF);
    sBI  = dmob * sBI;

    is_dir = false(size(sH)); lam = zeros(size(sH));

    if nblk == 1,
      [sF, sH, lam, is_dir] = handle_bc(sF, sH, lam, is_dir,   ...
                                        g, f_bc, h_bc, iF, iH, ...
                                        is_dirichlet,          ...
                                        sub_map.sub_f(face));
    end

    do_reg = ~any(is_dir);  % Need to set pressure zero level?
    [~, p, lam(~is_dir)] = schurComplementSymm(sBI, sC, sD(:,~is_dir), ...
                                               sF , sG, sH(  ~is_dir), ...
                                               'Regularize', do_reg, ...
                                               'LinSolve', linSolve);

    v = dmob * ([sC, sD] * [p; -lam]);   % v <- B*v == C*p - D*lam

    % Orthogonalize pressure against source term (theta(iG)) in each
    % coarse block.
    [p, press_ok(f)]= orth_press(p, theta(iG), cg.partition(iG), nBlk);

    % local numbering of cells and half-faces
    renum_p = zeros([niG,1],'uint32'); renum_p(:) = 1:niG;
    renum_v = zeros([niF,1],'uint32'); renum_v(:) = 1:niF;  

    % logical array. true for elements in iG that are in blk1
    in_blk1    = cg.partition(iG)==blk(1) ;
    % logical array. true for elements in iF that are in blk1
    hf_in_blk1 = ismembc( iF, sub_map.hf(iG(in_blk1)) );

    % new ordering
    ix_p = [ renum_p( in_blk1   ); ...
             renum_p(~in_blk1   ) ];    % cells from blk1 first
    ix_v = [ renum_v( hf_in_blk1); ...
             renum_v(~hf_in_blk1) ];    % hfs from blk1 first

    lpV(f,1) = {{iF(ix_v), v(ix_v), face, blk, ...
                 nhf_blk(blk), avgmob(blk), overlap}};
    lpP(f,1) = {{iG(ix_p), p(ix_p), face, blk, ...
                 ncells(blk) , avgmob(blk), overlap}};
  end

  % this will show only on workers where ~all(press_ok) == true
  % (avoids synching, of anything but std out ...)
  if ~all(press_ok)
    warning(['At least one pressure basis function does not ', ...
             'satisfy orthogonality condition. Quality of ' , ...
             'solution may be reduced.']);
  end

  % as long as dist is complete no communication is neccessary for this.
  coV = codistributed.build(lpV, dist, 'noCommunication');
  coP = codistributed.build(lpP, dist, 'noCommunication');
end
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function theta = get_coarse_weighting(g, p, w, src) % from MRST
   % All synthetic weighting terms must be strictly non-negative.
   %
   assert (~any(w < 0));

   % Initially, assume there are no explicit/external sources.
   %
   theta = w .* g.cells.volumes;

   if ~isempty(src),
      % Update for explicit sources if there nevertheless are some...

      % Determine coarse blocks already containing external sources.
      %
      has_src = accumarray(p(src.cell), 1, [max(p), 1]) > 0;

      % Eliminate previous (synthetic) weighting and apply correct source.
      %
      theta(has_src(p)) = 0;
      theta(src.cell)   = src.rate;
   end

   % Note:
   %   We need to normalize the (synthetic or explicit) source term 'theta'
   %   such that \int_{B_i} theta d\Omega == 1 lest the basis functions be
   %   inconsistent.
   %
   denom = accumarray(p, theta);  assert (all(abs(denom) > 0));
   theta = theta ./ denom(p);
end

%--------------------------------------------------------------------------

function [f, h, d] = expand_bc(g, bc) % adapted from MRST
   f = zeros([g.faces.num, 1]);
   h = zeros([g.faces.num, 1]);
   d = false([g.faces.num, 1]);

   if ~isempty(bc),
      assert (all(accumarray(bc.face, 1, [g.faces.num, 1]) <= 1));

      is_dir = strcmp('pressure', bc.type);
      f(bc.face(is_dir)) = bc.value(is_dir);
      d(bc.face(is_dir)) = true;

      is_neu = strcmp('flux', bc.type);
      h(bc.face(is_neu)) = - bc.value(is_neu);
   end
end

%--------------------------------------------------------------------------

function sub_c = sub_cells(g, p, overlap) % adapted from MRST
   nc    = g.cells.num;
   sub_c = sparse(1 : nc, p, 1, nc, max(p));  % == cg.cells.subCells

   if overlap > 0,
      n = double(g.faces.neighbors(all(g.faces.neighbors > 0, 2), :));
      n = sparse([n(:,1); n(:,2); (1 : nc).'], ...
                 [n(:,2); n(:,1); (1 : nc).'], 1, nc, nc);

      % BFS to discover immediate neighbours in overlap region.
      for o = 1 : overlap, sub_c = n * sub_c; end
   end

   sub_c = logical(sub_c);
end

%--------------------------------------------------------------------------

function mob = average_mobility(p, mob, vol)
   assert (all(numel(p) == [numel(mob), numel(vol)]));

   mob = accumarray(p, mob .* vol) ./ accumarray(p, vol);
end

%--------------------------------------------------------------------------

function [p, ok] = orth_press(p, w, b, nb) % from MRST
   present        = false([nb, 1]);
   present(b)     = true;

   renum          = zeros([nb, 1]);
   renum(present) = 1 : sum(present);

   % Compute orthogonality adjustment constant.  One value for each coarse
   % block present in the support of 'p'.
   a = accumarray(renum(b), p .* w) ./ ...
       accumarray(renum(b),      w);

   p = p - a(renum(b));

   % Check orthogonality to within relative bounds.
   ok = abs(p' * w) < 2 * numel(p) * eps(norm(p,inf));
end

%--------------------------------------------------------------------------

function [ops, nc, nhf] = mappings(g,cg,mob,bc,overlap) % from MRST
% mappings - Build essential grid mappings for coarse grid.

   % Compute 'sub_f' mapping depending on existence of any outer coarse
   % faces supporting flow (i.e., BC is not no-flow).  See 'help subFaces'
   % for details on MCOLON expression.
   if ~isempty(bc),
      [nsub, sub] = subFaces(g, cg);
      sub_ix      = cumsum([0; nsub]);
      sub_f       = @(f) sub(mcolon(sub_ix(f) + 1, sub_ix(f + 1)));
   else
      % Dangerous.  Assumes we won't be called upon to build a basis
      % function for an outer coarse face (f) if there are no external
      % boundary conditions...
      sub_f       = @(f) [];
   end

   % Count number of cells and number of half-faces in all coarse blocks.
   nc  = accumarray(cg.partition, 1);
   nhf = accumarray(cg.partition, double(diff(g.cells.facePos)));

   % sub_c(:,b) == true for all cells (rows) within (extended) block 'b'.
   sub_c = sub_cells(g, cg.partition, overlap);

   % cellno(i) == (fine-scale) cell which contains half-face 'i'.
   %
   cellno = rldecode(1 : g.cells.num, diff(g.cells.facePos), 2) .';

   %-----------------------------------------------------------------------
   % Define mapping operators. --------------------------------------------
   %
   %   1) ops.cells:
   %      any(sub_c(:,b), 2) is true for all fine-scale cells in the
   %      (possibly extended) block(s) 'b'.
   %
   ops.cells = @(b) find(any(sub_c(:,b), 2));

   %   2) ops.hf:
   %      Transpose because MCOLON gives row-vector, while we need columns.
   %
   ops.hf    = @(c) mcolon(g.cells.facePos(  c  ), ...
                           g.cells.facePos(c + 1) - 1) .';

   %   3) ops.faces:
   %      FIND those faces mentioned at least once in the half-faces 'i'.
   %
   ops.faces = @(i) find(accumarray(g.cells.faces(i,1), 1) > 0);

   %   4) ops.sub_f:
   %      See above.
   %
   ops.sub_f = sub_f;

   %   5) ops.dmob:
   %      mob(cellno(i)) is (total) mobility in cell containing half-face
   %      'i'.
   %
   ops.dmob  = @(i,n) spdiags(mob(cellno(i)), 0, n, n);
end

function [sF, sH, lam, is_dir] = handle_bc(sF, sH, lam, is_dir, ...
                                   g, f_bc, h_bc, iF, iH, ...
                                   is_dirichlet, ih) % from MRST

   fno     = zeros([g.faces.num, 1]);
   fno(iH) = 1 : numel(iH);

   is_bdry          = false(size(sH));
   is_bdry(fno(ih)) = true;
   is_dir(:)        = is_bdry & is_dirichlet(iH);
   is_neu           = is_bdry & ~is_dir;

   % This code assumes that a coarse face is either Dirichlet or
   % Neumann (not both).  That assumption must be revisited if the
   % following assertion fails.  We also fail if the face is neither
   % Dirichlet nor Neumann.
   %
   tot_flx = norm(h_bc(iH(is_neu)));
   assert (xor(sum(is_dir) > 0 && tot_flx < 1000 * eps(1), ...
               ~(tot_flx < 1000 * eps(1)) && sum(is_dir) == 0));

   % Dirichlet (pressure) boundary conditions.
   lam(is_dir) = f_bc(iH(is_dir));
   loc_iF      = is_dir(fno(g.cells.faces(iF,1)));
   sF(loc_iF)  = sF(loc_iF) - lam(is_dir);

   % Neumann (flux) boundary conditions.
   sH(is_neu)  = h_bc(iH(is_neu));
   denom = sum(sH(is_neu));
   if abs(denom) > sqrt(eps(denom)),
      sH(is_neu) = sH(is_neu) ./ denom;
   end
end
