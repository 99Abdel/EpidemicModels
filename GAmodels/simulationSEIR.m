clear all


N = 6e5;   % total number of population (if 1 is already adimensional)

I0 = 700;     % infected people (in the beginning) (if less than one is a portion of population)

a = 1.4;    % S to E coefficient (days^-1) (exposition rate)
b = .805;    % E to I coefficient (days^-1) (infection rate)
c = 0.61;   % I to R coefficient (days^-1) (removal rate)
