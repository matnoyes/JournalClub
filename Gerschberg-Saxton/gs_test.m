%% Test of GS from Thimons and Wittle
clear all; clc

% force a solution (1), or return first attempt result (0)?
force = 1;
% Track and animation convergence?
tr = 0;
% define sampling rate
delta = 1/256;
% max number of iterations
maxiter = 5000;
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


[estimate,track] = gs_tracking(f,F,maxiter,force);

fprintf('Number of Iterations: %d \n',estimate(3,1))
if estimate(2,length(estimate(2,:))) > 1e-3
    fprintf('Algorithm failed\n');
else
    fprintf('Solution Found\n');
    fprintf('Number of attempts: %d\n', estimate(4,1))
end

figure(7);
subplot(1,3,1); plot(t,phases); title('Original Phase'); ylabel('Phase [rad]');
subplot(1,3,2); plot(t,estimate(1,1:N)); title('Solution'); ylabel('Phase [rad]');
subplot(1,3,3); plot(1:estimate(3,1),estimate(2,1:estimate(3,1))); title('Error vs Iteration'); ylabel('|f''| - |g|');



if(tr)
    if estimate(3,1) == maxiter
        for j = 1:275
            if j < 6
            figure(8); plot(track(j,:)); pause(1);   title(['Iteration: ',num2str(j),' (',num2str(j/estimate(3,1)*100),'%)']);
            else
            figure(8); plot(track(j,:)); pause(.01); title(['Iteration: ',num2str(j),' (',num2str(j/estimate(3,1)*100),'%)']);
            end
        end
    else
        for j = 1:estimate(3,1)
            if j < 6
            figure(8); plot(track(j,:)); pause(1); title(['Iteration: ',num2str(j),' (',num2str(j/estimate(3,1)*100),'%)']);
            else
            figure(8); plot(track(j,:)); pause(.01); title(['Iteration: ',num2str(j)]);
            end
        end
    end
end

