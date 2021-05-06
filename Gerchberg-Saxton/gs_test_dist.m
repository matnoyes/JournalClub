%% A model of Gerchberg-Saxton's Probability Distribution Function for a given configuration
%  Matthew Noyes - JPL 2021
%  HCIT Journal Club
clear, close all; clc

NS = 10000; % Number of successed to run for the final Histogram
attempts = zeros(1,NS);

for l = 1:NS
    

    % define sampling rate
    samp = 1/256;
    % max number of iterations
    maxit = 500;
    % sample domain points in [ -.5 , .5)
    t = (-.5):samp:(.5-samp) ;
    % number of points ( used in FT calculations )
    N = length(t);

    % define rect function to limit domain
    syms x;
    r = piecewise(x<-.5, 0, x >= -.5 & x <=.5, 1, x >.5, 0);
    % sample object domain points
    full_f = double(subs(r,2*t)).*exp(30*1i*pi*t.^2);

    % calculate known phases for comparison
    phases = angle(full_f);
    % | f |
    f = abs(full_f);
    % | F |
    F = abs(fft(full_f)/N);


    [estimate,animate] = gs(f,F,maxit,1);

    p = l/NS*100;
    fprintf('NS Index: %d\nProgress: %4.2f%%\n\n',l,p)

    attempts(l) = estimate(4,1);

end

figure(8)
histogram(attempts);
title('# of failures before success (10,000 successes)');
xlabel('Number of sequential failures before success');
ylabel('count')