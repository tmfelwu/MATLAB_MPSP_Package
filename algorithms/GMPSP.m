function [t, X_opt, U_opt] = GMPSP(f,h, dF_dX, dF_dU, dH_dX, t0, tf, dt, X0, Yf, R, U_guess, tol, maxIter)
%GMPSP The function solves the optimal control algorithm.
%   f               : Dynamics of the system
%   h               : Output functoin y = H(x)
%   dF_dX           : Jacobian dF/dX
%   dF_dU           : Jacobian dF/dU
%   dH_dX           : Jacobian dH/dX
%   t0              : Initial Time
%   tf              : Final Time
%   dt              : Step Size
%   X0              : Initial State
%   Yf              : Final State
%   R               : Weigt matrix
%   U_guess         : Initial Guess with dimensions (m x N)
    N = floor(abs(tf - t0)/dt) + 1;
    U = U_guess;
    
    m = 1; 
    n = 2;
    p = 2;
    X = zeros(n,N);
    t = t0:dt:tf;
    
    X(:,1) = X0;
    
    for k_itr= 1:maxIter
        % Propagate System dynamics
        % ==============================
        for k = 1:(N-1)
            k1 = feval(f, t(k), X(:,k), U(:,k) );
            k2 = feval(f, t(k)+dt/2, X(:,k)+(dt/2)*k1, U(:,k) );
            k3 = feval(f, t(k)+dt/2, X(:,k)+(dt/2)*k2, U(:,k) );
            k4 = feval(f, t(k)+dt, X(:,k)+dt*k3, U(:,k) );

            X(:,k+1) = X(:,k) + (dt/6)*( k1 + 2*k2 + 2*k3 + k4 );
        end
        
        % Compute Output and the error in output
        % ================================
            Y_N = h(X(:,N));
            dY_N = Y_N - Yf;
            
        % Check for convergence
        % ================================
        if (k_itr > 1) % Allows at least one iteration

            if norm(Yf) > 1e-6
                err_Y_N = norm(dY_N)/norm(Yf);  % Percentage error
            else
                err_Y_N = abs(dY_N);     % Absolute error
            end
            if (err_Y_N < tol)
                disp('The algorithm has converged');
                k_itr;
                X_opt = X;
                U_opt = U;
                break;
            end
        end
        
        % Calculate the weight matrix W(t); dimension p x n 
        % =================================
        W = zeros(p,n,N);
        W_tf = dH_dX(X(:,end), U(:,end));
        W(:,:,end) = W_tf;
        for k = N:-1:2
            k1 = W(:,:,k) * dF_dX(X(:,k), U(:,k));
            k2 = (W(:,:,k) + (dt/2) * k1 ) * dF_dX(X(:,k), U(:,k));
            k3 = (W(:,:,k) + (dt/2) * k2 ) * dF_dX(X(:,k), U(:,k));
            k4 = (W(:,:,k) + dt * k3 ) * dF_dX(X(:,k), U(:,k));
            W(:,:,k-1) = W(:,:,k) + (dt/6)*( k1 + 2*k2 + 2*k3 + k4 );      
        end
        
        % Calculate B_s(t) p * m
        % =================================
        B_s = zeros(p,m,N);
        for k=1:N
            B_s(:,:,k) = W(:,:,k) * dF_dU(X(:,k), U(:,k));        
        end        
        % Calculate A_lambda and B_lambda
        % =================================
        int1 = zeros(p,p,N);
        for k=1:N
            int1(:,:,k) = B_s(:,:,k) * inv(R) * B_s(:,:,k)';
        end
        A_lambda =  trapz(t0:dt:tf, int1, 3);
        int2 = zeros(p,1,N);
        for k=1:N
            int2(:,:,k) = B_s(:,:,k) * U(:,k);
        end
        B_lambda =  trapz(t0:dt:tf, int2, 3);
        % Calculate the new Control
        % =================================
        for k =1:N
            U(:,k) = - inv(R) * B_s(:,:,k)' *  inv(A_lambda) * ( dY_N  - B_lambda);
        end
        
        if (k_itr >= maxIter)
            disp('Maximum number of iterations reached.');
            disp('The algorithm has not converged');
            X_opt = X;
            U_opt = U;
        end
    end
end

