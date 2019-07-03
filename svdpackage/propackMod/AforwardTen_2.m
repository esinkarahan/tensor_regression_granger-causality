function y = AforwardTen_2(x)
% Auxuliary function used for PROPACK LANEIG
global AA
At = AA';
y = At*x;
