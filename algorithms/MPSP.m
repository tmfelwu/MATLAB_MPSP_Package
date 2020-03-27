function [t_opt, U_opt, X_opt, Y_N, dY_N, k_itr] = MPSP(f, h, Jf_X, Jf_U, Jh_X, t0, tf, dt, X0, Y_N_st, R, U_guess, tol_Y_N, max_itr, option)
%
% This function, developed by Prof. Radhakant Padhi, implements the MPSP 
% algorithm to find optimal control U_opt and optimal state X_opt trajectories.
% It also outputs the actual Y_N as well as the error of output
% from the desired value Y_N_st, i.e. dY_N (in absolute terms). In
% addition, it also outputs the number of iterations it takes to converge, i.e. k_itr.
%
% Usage: [U_opt, X_opt, Y_N, dY_N, k_itr] = MPSP(f, h, Jf_X, Jf_U, Jh_X, t0, tf, dt, X0, Y_N_st, R, U_guess, tol_Y_N, max_itr, option)
%
% The necessary inputs are:
% (1)  The state function  'f' in X_dot = f(X,U)
% (2)  The output function 'h' in Y = h(X)
% (3)  The Jacobian matrix 'Jf_X' = df/dX
% (4)  The Jacobian matrix 'Jf_U' = df/dU
% (5)  The Jacobian matrix 'Jh_X' = dh/dX
% (6)  Initial time t0
% (7)  Final time tf
% (8)  Sampling time dt
% (9)  Initial condition X0
% (10) Desired final value of the output vector Y_N_st
% (11) Weighting matrix R (must either be m x m constant matrix or a m x m x (N-1) 3D matrix)
% (12) Control guess history U_guess [dimension should be m x (N-1), where N is the
%      number of time grids]
% (13) Tolerence value for termination (optional). Default value = 2% (in 2-Norm sense)
% (14) Maximum number of iterations in case the algorithm does not converge
%      (optional). Default value = 50
% (15) Option = 1 (control mimimization, default), Option = 2 (control
%      deviation mimimization (from previous iteration)
%
% Copyright: Prof. Radhakant Padhi, Indian Institute of Science, Bangalore,
% India. All rights reserved. This function is not for distribution without
% prior notice and approval from Prof. Radhakant Padhi.

if (nargin == 12)
    tol_Y_N = 0.02;  % Default value for output tolerance. This is
    % usually in percentage error sense, if Y_N_st is not
    % close to zero. Else in absolute error sense (see for
    % details later in the code)
    max_itr = 50;    % Default value for number of maximum number of iterations
    option = 1;      % Defaul for control mminimization
end

N = floor(abs(tf - t0)/dt) + 1;

m_X = length(X0);  % m_X is required later
[m_U, N_U] = size(U_guess);
[m_R, n_R, N_R] = size(R);

if (N_U ~= (N-1))
    error('Dimension of control history guess is not compatible with the time grid');
end

if ((m_R ~= n_R) | (m_R ~= m_U))
    error('Dimension of control weightage matrix is not compatible');
end

if ((N_R ~= 1) & (N_R ~= (N-1)))
    error('Dimension of control weightage matrix is not compatible');
end

if (N_R == 1)
    for k = 1:(N-1)
        R_k(:,:,k) = R;   % Constant weightage
    end
elseif (N_R == (N-1))
    R_k = R;
end

% Initialization1
% ==============
t(:,1) = t0;
X(:,1) = X0;
U = U_guess;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for (k_itr = 1:max_itr)

    % State propagation (using RK_4 with fixed step size)
    % ===================================================
    for k = 1:(N-1)

        k1 = feval(f, t(k), X(:,k), U(:,k) );
        k2 = feval(f, t(k)+dt/2, X(:,k)+(dt/2)*k1, U(:,k) );
        k3 = feval(f, t(k)+dt/2, X(:,k)+(dt/2)*k2, U(:,k) );
        k4 = feval(f, t(k)+dt, X(:,k)+dt*k3, U(:,k) );
        % Note: U(:,k) comes as a parameter since it remains constant within
        % the interval between k and (k+1)

        X(:,k+1) = X(:,k) + (dt/6)*( k1 + 2*k2 + 2*k3 + k4 );

        t(:,k+1) = t(:,k) + dt;
    end

    % Evaluation of output and output error
    % =====================================

    Y_N = h(X(:,N));

    dY_N = Y_N - Y_N_st;

    % Termination of the algorithm
    % =============================

    if (k_itr > 1) % Allows at least one iteration

        if norm(Y_N_st) > 1e-6
            err_Y_N = norm(dY_N)/norm(Y_N_st);  % Percentage error
        else
            %             err_Y_N = dY_N*100;     % Absolute error
            %             % Amplification of actual error (for termination condition)
            %             % ensures lesser final error
            err_Y_N = abs(dY_N);     % Absolute error
        end


        if (err_Y_N < tol_Y_N)
            disp('The algorithm has converged');
            k_itr
            X_opt = X;
            U_opt = U;
            break;
        end
    end

    % Evaluation of Sensitivity Matrices
    % ==================================

    for (k = 1:(N-1))
        df_dX(:,:,k) = feval(Jf_X, X(:,k), U(:,k));
        dF_dX(:,:,k) = eye(m_X) + df_dX(:,:,k) * dt; % With Euler discretization

        df_dU(:,:,k) = feval(Jf_U, X(:,k), U(:,k));
        dF_dU(:,:,k) = df_dU(:,:,k) * dt;            % With Euler discretization
    end

    dY_N_dX_N = feval(Jh_X, X(:,N));    % This is (dY_N/dX_N)
    B_0(:,:,N-1) = dY_N_dX_N;

    for (k = (N-2) : -1 : 1)
        B_0(:,:,k) = B_0(:,:,k+1) * dF_dX(:,:,k+1);
    end

    for (k = (N-1): -1 : 1)
        B(:,:,k) = B_0(:,:,k) * dF_dU(:,:,k);
    end

    % Control history update
    % ======================

    p = length(Y_N);
    A_lambda = zeros(p,p);
    for k = 1:(N-1)
        A_lambda = A_lambda - B(:,:,k)*inv(R_k(:,:,k))*B(:,:,k)';
    end

    if (option == 1) % Control Minimization

        b_lambda = zeros(p,1);
        for k = 1:(N-1)
            b_lambda = b_lambda + B(:,:,k)*U(:,k);
        end

        lambda = -A_lambda \ (dY_N - b_lambda);

        for k = 1:(N-1)
            dU(:,k) = R_k(:,:,k) \ (B(:,:,k)'*lambda) + U(:,k);
        end

    elseif (option == 2)    % Control Deviation Minimization

        lambda = A_lambda \ dY_N;

        for k = 1:(N-1)
            dU(:,k) = -R_k(:,:,k) \ (B(:,:,k)'*lambda);
        end

    end

    U_new(:,:,k_itr) = U - dU; % To see how the control history update proceeds

    U = U_new(:,:,k_itr);
    t_opt = t(1:N_U);

    if (k_itr >= max_itr)
        disp('Maximum number of iterations reached.');
        disp('The algorithm has not converged');
        X_opt = X;
        U_opt = U;
    end
end
