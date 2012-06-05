function [g,rock,cg,mob,bc,...
         src,overlap,facetrans,...
         weighting,activeBnd] = ...
         xBroadcast(g,rock,cg,mob,bc,...
         src,overlap,facetrans,...
         weighting,activeBnd)

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