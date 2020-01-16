% Given a water height and a z coordinate system,
% which nodes have water and which do not?
%   wet = 1: node touching moulin water
%   wet = 0: node touching air
function wet = locatewater(hw,z)
%
wet = z<=hw;