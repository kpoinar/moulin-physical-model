% Given a water height and a z coordinate system,
% which nodes have water and which do not?
%   wet = 1: node touching moulin water (wet)
%   wet = 0: node touching air          (not wet)
function wet = locatewater(hw,z)
%
wet = z<=hw;