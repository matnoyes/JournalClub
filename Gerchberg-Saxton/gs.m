function estimate = gs(f,F,maxiter,force)
% param f : | sampled points |
% param F : | FT of sampled points |
% param phi : intial phase estimate
% param maxiter : max number of iterations
% param force: Will determine if the function will re-run until a solution
% is found, or just return the failed attempt
% return estimate : 3 row array 
% (phase estimates, iterativeerror, # of iterations tosuccess or max)


if force == 1
    
solution = 0;
attempts = 0;

while (solution == 0)

% random phases to start algorithm
N = length(f);
phi = 2*pi*rand(1, N)-pi;
% number of points ( used in FT and error calculations)
N = length(f);
% define error tolerance for success
error_tol = 1e-3;
% to keep track of error through iterations
error = zeros(2,maxiter);
% 1 st estimate with given phases
x = f.*exp(1i.*phi);

k =1;
% repeat algorithm until convergence or max number of iterations
while (k == 1 || (error(k-1) > error_tol && k < maxiter +1))
% eq 5
X = fft( x )/N ;
% Fourier domain error
error(1,k) = sqrt(N)*norm(abs(X) - F);
% eq 6
Y = F.*exp(1i.*angle(X));

% eq 7
y = N*ifft(Y);
% object domain error
error(2,k) = norm(abs(y) - f );
% eq 4
x = f.*exp(1i.*angle(y));

% increment index
k = k + 1;


end

if error(k-1) > error_tol %star over
    clear N error x k X error Y y phi;
elseif (error(k-1) <= error_tol)
    solution = 1;
end

attempts = attempts+1;
end

% ending estimated phases
est_phases = angle(x);
estimate(1, 1: length(est_phases)) = est_phases;
% combine errors into single array
estimate(2, 1:2*(k-1)) = reshape(error(1:2, 1:k -1), 1, 2*(k-1));
% number of iterations it took to reach success / max
estimate(3,1) = k-1;
estimate(4,1) = attempts;

elseif force == 0
    % random phases to start algorithm
N = length(f);
phi = 2*pi*rand(1, N)-pi;
% number of points ( used in FT and error calculations)
N = length(f);
% define error tolerance for success
error_tol = 1e-3;
% to keep track of error through iterations
error = zeros(2,maxiter);
% 1 st estimate with given phases
x = f.*exp(1i.*phi);

k =1;
% repeat algorithm until convergence or max number of iterations
while (k == 1 || (error(k-1) > error_tol && k < maxiter +1))
% eq 5
X = fft( x )/N ;
% Fourier domain error
error(1,k) = sqrt(N)*norm(abs(X) - F);
% eq 6
Y = F.*exp(1i.*angle(X));

% eq 7
y = N*ifft(Y);
% object domain error
error(2,k) = norm(abs(y) - f );
% eq 4
x = f.*exp(1i.*angle(y));

% increment index
k = k + 1;

% ending estimated phases
est_phases = angle(x);
estimate(1, 1: length(est_phases)) = est_phases;
% combine errors into single array
estimate(2, 1:2*(k-1)) = reshape(error(1:2, 1:k -1), 1, 2*(k-1));
% number of iterations it took to reach success / max
estimate(3,1) = k-1;
estimate(4,1) = 1;
end
end