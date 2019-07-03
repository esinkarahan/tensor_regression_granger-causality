function y = Aforward_2(x)
% Auxuliary function used for PROPACK LANEIG
global A
At = A';
y = At*x;
