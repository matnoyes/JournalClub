%% Test of GS from Thimons and Wittle
clear, close all; clc

for l = 1:10000
 
%     i
%     clear estimate
% define sampling rate
delta = 1/256;
% max number of iterations
maxiter = 500;
% sample domain points in [ -.5 , .5)
t = (-.5):delta:(.5-delta) ;
% number of points ( used in FT calculations )
N = length(t);

% define rect function to limit domain
syms x;
r = piecewise(x<-.5, 0, x >= -.5 & x <=.5, 1, x >.5, 0);
% sample object domain points
full_f = double(subs(r,2*t)).*exp(30*i*pi*t.^2);

% calculate known phases for comparison
phases = angle(full_f);
% | f |
f = abs(full_f);
% | F |
F = abs(fft(full_f)/N);



% figure(1);
% plot(t,exp(30*i*pi*t.^2))
% 
% figure(2);
% plot(t,full_f)
% 
% figure(3);
% plot(t,phases)
% 
% figure(4);
% plot(t,f);
% 
% figure(5);
% plot(fftshift(F))
%
% figure(6);
% plot(phi)


estimate = gs(f,F,maxiter);

% 
% fprintf('Number of Iterations: %d \n',estimate(3,1))
% if estimate(2,length(estimate(2,:))) > 1e-3
%     fprintf('Algorithm failed\n');
% else
%     fprintf('Solution Found\n');
%     fprintf('Number of attempts: %d\n', estimate(4,1))
% end

% figure(7);
% subplot(1,3,1); plot(1:estimate(3,1),estimate(2,1:estimate(3,1))); title('Error vs Iteration');
% subplot(1,3,2); plot(t,estimate(1,1:N)); title('Original and Solution');
% subplot(1,3,3); plot(t,phases); title('Original Phase');

l
attempts(l) = estimate(4,1);
end

figure(8)
histogram(attempts);
title('# of failures before success (10,000 successes)');
xlabel('Number of sequential failures before success');
ylabel('count')