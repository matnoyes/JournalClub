function [estimate,animate] = gs(f,F,maxiter,force)
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
        % Random phases to start algorithm
        N = length(f);
        phi = 2*pi*rand(1,N)-pi;
        
        N = length(f); % number of points (used in FT and error calculations)
        
        error_tol = 1e-3; % define error tolerance for success
        
        error = zeros(2,maxiter); % to keep track of error through iterations
        
        x = f.*exp(1i.*phi); % 1 st estimate with given phases

        animate(1,1:N) = angle(x); % 1st estimate recorded
        
        k =1;
        
        % repeat algorithm until convergence or max number of iterations
    while (k == 1 || (error(k-1) > error_tol && k < maxiter +1))
        
        X = fft( x )/N ;
        
        error(1,k) = sqrt(N)*norm(abs(X) - F); % Fourier domain error
        
        Y = F.*exp(1i.*angle(X));

        y = N*ifft(Y);
        
        error(2,k) = norm(abs(y)-f); % object domain error
        
        x = f.*exp(1i.*angle(y));

        k = k + 1;
        
        animate(k,1:N) = angle(x); % Subsequent estimates recorded

    end

    if error(k-1) > error_tol % Start over if no solution found
        clear N error x k X error Y y phi;
    elseif (error(k-1) <= error_tol)% End if solution found
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
    % number of attempts it took to find solution
    estimate(4,1) = attempts;

elseif force == 0 % This block will run the algorithm once, and will return an answer regardless of whether or not its a solution
    % Random phases to start algorithm
        N = length(f);
        phi = 2*pi*rand(1,N)-pi;
        
        N = length(f); % number of points (used in FT and error calculations)
        
        error_tol = 1e-3; % define error tolerance for success
        
        error = zeros(2,maxiter); % to keep track of error through iterations
        
        x = f.*exp(1i.*phi); % 1 st estimate with given phases

        animate(1,1:N) = angle(x); % 1st estimate recorded
        
        k =1;
        
        % repeat algorithm until convergence or max number of iterations
    while (k == 1 || (error(k-1) > error_tol && k < maxiter +1))

        X = fft( x )/N ;
        
        error(1,k) = sqrt(N)*norm(abs(X) - F); % Fourier domain error
        
        Y = F.*exp(1i.*angle(X));

        y = N*ifft(Y);
        
        error(2,k) = norm(abs(y)-f); % object domain error
        
        x = f.*exp(1i.*angle(y));

        k = k + 1;
        
        animate(k,1:N) = angle(x); % Subsequent estimates recorded

        % ending estimated phases
        est_phases = angle(x);
        estimate(1, 1: length(est_phases)) = est_phases;
        % combine errors into single array
        estimate(2, 1:2*(k-1)) = reshape(error(1:2, 1:k -1), 1, 2*(k-1));
        % number of iterations it took to reach success / max
        estimate(3,1) = k-1;
        % number of attempts it took to find solution
        estimate(4,1) = 1;
    end
end