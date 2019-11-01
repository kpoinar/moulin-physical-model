function [h,S,Qout] = subglacialsc(Mr,R,)

% This function is based on Schoof 2010, without the cavity.
% input:
% Mr = Moulin radius in function of z. 
% R = supraglacial recharge. Can be constant or a function of time.
% 
% output:

% import constants
C = makeConstants;
