function [g,rock,cg,mob,bc,src,overlap,facetrans,...
          weighting,activeBnd] = xInitWorkers(p)

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
