%% Implementation of the Stochastic Smoothing Phase Retrieval algorithm

function [z0,z,Relerrs] = SSPR(x,y, Params, A, At,Amatrix)    
%% Initialization
    
    z0      = randn(Params.n, 1);
    z      = z0 / norm(z0, 'fro');                          % Initial guess
    normest = sqrt(sum(y(:)) / numel(y(:)));                % Estimate norm to scale eigenvector
    m       = Params.m;
    ymag    = sqrt(y);

    [ysort, ~] = sort(y, 'ascend');
    ythresh    = ysort(round(m / 1.3));
    ind        = y >= ythresh;

    Aselect    = Amatrix(ind, :);
    weights    = (ymag(ind)).^(Params.alpha);               % weights w_i

    %% The weighted maximal correlation initialization can be computed using power iterations
    %% or the Lanczos algorithm, and the latter performs well when m/n is small
    for tt = 1:Params.npower_iter                           % Power iterations
        z  = Aselect' * (weights .* (Aselect * z));
        z  = z / norm(z, 'fro');
    end

    z = normest * z;
    Relerrs = norm(x - exp(-1i * angle(trace(x' * z))) * z, 'fro') / norm(x, 'fro'); % Initial rel. error
    
    u = Params.u0;
    %% main loop
    for t = 1: Params.T*Params.m,
        tt = mod(t,Params.m)+1;
        grad = compute_grad(z,u,ymag,Amatrix,tt);
        
        z = z-Params.mu*grad;
        
        if norm(compute_grad(z,u,ymag,Amatrix,tt),'fro')< Params.y*u
            u = Params.y1*u;
        end
        
        if tt==1
            Relerrs = [Relerrs, norm(x - exp(-1i*angle(trace(x'*z))) * z, 'fro')/norm(x,'fro')]; % Relative error
        end
    end
