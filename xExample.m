% code assumes MRST is already initiated

% helper function to initiate workers 
%   and construct composite variables
% all returned variables are Composite
% the following variables are assigned default values:
%   overlap   = 0;
%   bc        = [];
%   src       = [];
%   facetrans = zeros(0,2);
%   activeBnd = [];
%   weighting = 'perm';
% the remaining variables are empty.
% they need to be initialized as in the example below. 
% afterwards xBroadcast needs to be called.
procs = 4 ; % 4 workers
[g,rock,cg,mob,bc,src,overlap,facetrans,weighting,...
 activeBnd] = xInitWorkers(procs);

% disable gravity on all workers
pctRunOnAll gravity off;

% parallel block. 
% variables that define the problem are initialized on
% worker 1
% they are then broadcasted using xBroadcast.
spmd
  % initialization of the problem only on worker 1
  if labindex == 1,

    fprintf('initializing on worker 1 ... ');

    % size of fine and coarse geometry
    nx = 100; ny = 50; nz = 15;
    Nx = 10;  Ny = 5;  Nz = 5;
    
    % geometry and rock properties
    g = computeGeometry(processGRDECL(makeModel3(...
      [nx, ny, nz])));
    K = logNormLayers(g.cartDims, [10, 300, 40, 0.1, 100]);
    rock.perm = bsxfun(@times, [1, 100, 0.1], K(:));
    rock.perm = convertFrom(rock.perm(g.cells.indexMap, :),...
      milli*darcy);

    % coarse partition
    p  = partitionUI(g, [Nx, Ny, Nz]);
    p  = processPartition(g,p);
    cg = generateCoarseGrid(g,p);
    
    % fluid properties
    fluid = initSingleFluid('mu' ,1*centi*poise, ...
                            'rho',1000*kilogram/meter^3);

    % solution structure
    xMs = initState(g, [], 0, [0, 1]);

    % get mobility
    mu  = fluid.properties(xMs);
    kr  = fluid.relperm(ones([g.cells.num,1]),xMs);
    mob = kr ./ mu;
    
    % some driving forces.
    % sources can be added in a similar fashion here
    bc = pside(bc,g,'East',1);
    bc = pside(bc,g,'West',0);

    fprintf('done\n');
    
    % should you wish to call functions 
    %   from regular MRST for comparison
    %   they can be called on one worker 
    %   inside spmd blocks such as this,
    %   using the same composite variables; 
    %   eg:
    %   
    %   ss = computeMimeticIP(g,rock);
    %   cs = generateCoarseSystem(g,rock,ss,cg,mob,'bc',bc);    
  end

  % distribute Composite variables to *all* workers
  if labindex == 1, tic; end % time on worker 1.
  [g,rock,cg,mob,bc,src,overlap,facetrans,...
   weighting,activeBnd] = xBroadcast(...
   g,rock,cg,mob,bc,src,overlap,facetrans,...
   weighting,activeBnd);
  if labindex == 1,
    fprintf('time xBroadcast:\t%.5f\n',toc);
  end
end

% calculate inner products.
% - XBI -- distributed cell array. 
%          XBI{i} contains inner product of cell i.
tic;
XBI = xComputeMimeticIP(g,rock,facetrans);
fprintf('time xComputeMimeticIP:\t%.5f\n',toc);

% distribution of inner products
% - coX     -- distributed cell array.
%                contains inner products
% - coiG    -- distributed cell array.
%                contains global numbering of cells
% - coFaces -- distributed array. 
%                contains global numbering of faces
tic;
[coXBI,coiG,coFaces] = xDistributeIP(XBI,g,cg,overlap,...
                                     bc,activeBnd);
fprintf('time xDistributeIP:\t%.5f\n',toc);

% construct basis functions
% - xCS -- composite array
%            non-empty only on worker 1. 
%            contains coarse basis functions.
tic;
XCS = xGenerateCoarseSystem(g,rock,cg,overlap,...
                            bc,src,weighting,mob,...
                            coXBI,coiG,coFaces);
fprintf('time xGenerateCoarseSystem:\t%.5f\n',toc);

% build BI and S from X, and solve coarse system:
spmd
  BI = gather(XBI,1);
  if labindex == 1
    dimProd = double(diff(g.cells.facePos));
    [ind1, ind2] = blockDiagIndex(dimProd, dimProd);
    n = size(g.cells.faces, 1);
    S.BI = sparse(ind1,ind2,vertcat(BI{:}),n,n);
    S.type = 'hybrid'; S.ip = 'ip_simple';
    
    xMs = solveIncompFlowMS(xMs,g,cg,p,S,XCS,fluid,...
      'bc',bc,'Solver',S.type);
  end
end

% can not plot from workers
% transport xMs and g to host
xxMs = xMs{1}; gg = g{1};
% note that composite variables do not 
% support accessing structure fields
% using '.', such as myStruct.field
clf,
  plotCellData(gg, convertTo(xxMs.pressure, barsa()));
  view(3), camproj perspective, axis equal tight off, 
  camlight headlight
  cax = caxis; cobar =  colorbar;
