%% A model of Gerchberg-Saxton
%  Matthew Noyes - JPL 2021
%  HCIT Journal Club
clear all; clc

% Force a solution (1), or return first attempted result regardless of success (0)?
force = 0;
% Track and animate convergence [on (1) / off (0)]?
ani = 0;


% Define sampling rate
samp = 1/256;
% Max number of iterations to run before failure is decided
maxit = 1000;
% Sample domain points in (-.5, .5)
t = (-.5):samp:(.5-samp) ;
% Number of points (used in FT calculations)
N = length(t);

% Define rect function to limit domain
syms x;
r = piecewise(x<-.5, 0, x >= -0.5 & x <= 0.5, 1, x >.5, 0);
% Sample object domain points
f_full = double(subs(r,2*t)).*exp(30*i*pi*t.^2);

% Calculate known phases for comparison
phases = angle(f_full);
% | f |
f = abs(f_full);
% | F |
F = abs(fft(f_full)/N);

[estimate,animate] = gs(f,F,maxit,force);

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


if(ani) % This block animates the convergence
    if estimate(3,1) == maxit % If a solution was NOT found, just animate the first 275 iterations
        for j = 1:275
            if j < 6 % Animate the first 6 iterations slowly
            figure(8); plot(animate(j,:)); pause(1);   title(['Iteration: ',num2str(j),' (',num2str(j/estimate(3,1)*100),'%)']);
            else     % Then speedrun the rest of the animation
            figure(8); plot(animate(j,:)); pause(.01); title(['Iteration: ',num2str(j),' (',num2str(j/estimate(3,1)*100),'%)']);
            end
        end
    else
        for j = 1:estimate(3,1) %If a solution WAS found, animate all iterations
            if j < 6 % Animate the first 6 iterations slowly
            figure(8); plot(animate(j,:)); pause(1); title(['Iteration: ',num2str(j),' (',num2str(j/estimate(3,1)*100),'%)']);
            else     % Then speedrun the rest of the animation
            figure(8); plot(animate(j,:)); pause(.01); title(['Iteration: ',num2str(j)]);
            end
        end
    end
end

