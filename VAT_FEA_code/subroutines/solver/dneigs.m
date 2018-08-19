function varargout = dneigs(varargin)
% dneigs find a few eigenvalues and eigenvectors of a given large sparse
% nonsymmetric square matrix A of order n. 
%
% dneigs solve the standard eigenvalue problem A*v = lambda*v. 
%
% d = dneigs(A) returns a vector of A's 6 largest magnitude eigenvalues.
%
% [v,d] = dneigs(A) returns a matrix v whose columns are the eigenvectors
% corresponding to the diagonal matrix d of A's 6 largest magnitude
% eigenvalues.
%
% [v,d,ds,flag] = dneigs(A) returns a convergence flag. If the value of
% flag is 0 then all the requested number of eigenvalues converged;
% otherwise not all converged. ds display the matrix-vector products,number
% of iteration,number of orthogonalization steps.
%
% In this version, we consider only the case when A is a given matrix.
% Later we will concentrate on the case when A is a function that returns
% the vector y = A*x.
%
% dneigs(A,k) returns the k largest magnitude eienvalues.
%
% dneigs(A,k,sigma) returns k eigenvalues with respect to sigma. We consider
% the cases when sigma is string only. Later we will consider the case when
% sigma is an scalar. Sigma takes the following forms:
%
% 'LM' or 'SM' - Largest or Smallest Magnitude.
% 'LR' or 'SR' - Largest or Smallest Real part.
%
% dneigs(A,k,sigma,opts) specifies options, which must be a structure:
%
% opts.tol: convergence tolerance: norm(A*v-v*d)<=tol [scalar | {1e-6}]
% opts.p_min: minimum number of Arnoldi vectors: k+1<p_min<=n [integer | {k+2}]
% opts.p_max: maximum number of Arnoldi vectors: p_min<p_max<=n [integer | {2*k+4}]
% opts.maxit: maximum number of iterations [integer | {300}]
% opts.disp: diagnostic information display level [0 | {1} | 2]
% opts.v0: starting vector [n-by-1 vector | {rand(n,1)-0.5}]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Main Routine                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% processing the inputs and do error-checking
if (nargin > 4),
    error('Too many input arguments.')
end

if (nargout > 4),
    error('Too many output arguments.')
end

% checking the first input argument
if isa(varargin{1}, 'double'),
    A = varargin{1};
    Amatrix = true;
else
    error('Matrix A must be double, neigs does not support AFUN.')
end

if Amatrix,
    [m,n] = size(A);
    if (m ~= n),
        error('Matrix A must be square.')
    end
end

isrealprob = 1; % isrealprob = isreal(A)
if Amatrix,
    isrealprob = isreal(A);
end

if ~isstr(A),
    pmmd = symamd(A);
    A = A(pmmd,pmmd);
end

% checking the second input argument
if (nargin < (4 - Amatrix - 1)),
    k = min(n,6); % setting default value of k
else 
    k = varargin{4 - Amatrix - 1};
end

kstr = ['The requested number of eigenvalues, k, must be a positive ' ...
    'integer and satisfies the constraint k <= n.'];

if (~isa(k,'double') || ~isscalar(k) || ~isreal(k) || (k > n)),
    error(kstr)
end

if issparse(k),
    k = full(k);
end

if (round(k) ~= k),
    warning('Rounding the number of eigenvalues.', kstr) 
end

% trick to get faster convergence, the tail end (closest to cut off of the 
% sort) will not converge as fast as the leading end of the "wanted" Ritz 
% values 
ksave = k;
% k = min(n,k + 6);
ksavesave = ksave;

% checking the third input argument
sigmastr = sprintf(['The nature of the eigenvalue, sigma, must be a ' ...
    'valid 2-element string.']);

if (nargin < (5 - Amatrix - 1)),
    sigma = 'LM'; % setting default value of sigma
else
    sigma = varargin{5 - Amatrix - 1};
    if (~ischar(sigma) || ~isequal(size(sigma),[1,2])),
        error(sigmastr, 'the choices are LM, SM, LR, or SR.')
    end
   
    sigma = upper(sigma);
    
    if ~ismember(sigma, {'LM', 'SM', 'LR', 'SR'}),
        error(sigmastr, 'the choices are LM, SM, LR, or SR.')
    end
end

% checking the fourth input argument
if (nargin >= (6 - Amatrix - 1)),
    opts = varargin{6 - Amatrix - 1};
    if ~isa(opts, 'struct'),
        error('The options argument must be of structure form.')
    end
else 
    opts = [];
end

% assiging the variables
if isfield(opts,'tol'),
    if (~isequal(size(opts.tol),[1,1]) || ~isreal(opts.tol) || ...
            (opts.tol <= 0)),
        error(['The tolerance for convergence must be a strictly ' ...
            'positive real scalar.'])
    else
        tol = full(opts.tol);
    end
else
    % setting default value for tol
    tol = 1e-6;
end

if isfield(opts,'p_min'),
    pstr = ['The minimum number of Arnoldi vectors opts.p_min must be a ' ...
        'positive integer <= n'];
    if (~isequal(size(opts.p_min),[1,1]) || ~isreal(opts.p_min) || ...
            (opts.p_min <= 0) || (opts.p_min > n)),
        error(pstr)
    else
        p_min = opts.p_min;
    end

    if issparse(p_min),
        p_min = full(p_min);
    end

    if (round(p_min) ~= p_min),
        warning('Rounding the value of p_min', pstr)
        p_min = round(p_min);
    end

    if (p_min <= k + 1),
        error('The minimum number of Arnoldi vectors opts.p_min > k + 1.')
    end
else
    % setting default value of p_min
    if (k < 8),
        p_min = k + 2;
    elseif (k < 20),
        p_min = k + 2;
    else
        p_min = k + 2;
    end
end

if isfield(opts,'p_max'),
    pstr = ['The maximum number of Arnoldi vectors opts.p_max must be a ' ...
        'positive integer <= n'];
    if (~isequal(size(opts.p_max),[1,1]) || ~isreal(opts.p_max) || ...
            (opts.p_max <= 0) || (opts.p_max > n)),
        error(pstr)
    else
        p_max = opts.p_max;
    end

    if issparse(p_max),
        p_max = full(p_max);
    end

    if (round(p_max) ~= p_max),
        warning('Rounding the value of p_max', pstr)
        p_max = round(p_max);
    end

    if (p_max <= p_min),
        error('The maximum number of Arnoldi vectors opts.p_max > p_min.')
    end
else
    % setting default value of p_max
    if (k < 8),
        p_max = 19;
    elseif (k < 30),
        p_max = 3*k - fix(k/2);
    else
        p_max = 3*k - fix(3*k/6);
    end
end

if isfield(opts,'maxit'),
    str = ['The required number of maximum iteration opts.maxit ' ...
        'must be a positive integer.'];
    if (~isequal(size(opts.maxit),[1,1]) || ~isreal(opts.maxit) || ...
            (opts.maxit <= 0)),
        error(str)
    else
        maxit = opts.maxit;
    end

    if issparse(maxit),
        maxit = full(maxit);
    end

    if (round(maxit) ~= maxit),
        warning('Rounding the value of maxit',str)
        maxit = round(maxit);
    end
else
    maxit = max(300,ceil(2*n/max(p_max,1)));
end

if isfield(opts,'disp'),
    dispstr = ['Diagnostic level opts.disp must be an integer.'];
    if (~isequal(size(opts.disp),[1,1]) || ~isreal(opts.disp) || ...
            (opts.disp < 0)),
        error(dispstr)
    else
        display = opts.disp;
    end

    if (round(display) ~= display),
        warning('Rounding the value disp',dispstr)
        display = round(display);
    end
else
    display = 1;
end

if isfield(opts, 'v0'),
    if ~isequal(size(opts.v0),[n,1])
        error('The starting vector v0 must be a vector of length n.')
    end

    if isrealprob,
        if ~isreal(opts.v0),
            error(['For reak problem the starting vector opts.v0 ' ...
                'must be real.'])
        end
        v0 = full(opts.v0);
    else
        v0 = full(opts.v0);
    end
else
    %setting default starting vector
    rand('seed',0)
    v0 = rand(n,1) - 0.5;
end

% start timing
t0 = cputime;

% assigning variables
info = 0;
v = v0;
matvec = 0;
ortcost = 0;
matvectime = 0;
orttime = 0;
beta = 1.0/norm(v(:,1));
v(:,1) = v(:,1)*beta;

% computing the matrix-vector product w = Av
v(:,2) = A*v(:,1);

% computing the first entry of the upper Hessenberg matrix H
alpha = v(:,1)'*v(:,2);
h(1,1) = alpha;

% computing the residual vector in the second column of the matrix V
v(:,2) = v(:,2) - alpha*v(:,1);

% performing one step of reorthogonalization to correct any
% orthogonality problems
alpha = v(:,1)'*v(:,2);
v(:,2) = v(:,2) - alpha*v(:,1);
h(1,1) = h(1,1) + alpha;

% Now we can compute k-steps of the Arnoldi sequence
% assigning the variables
kstart = 1;
ritz = 1.0;
kp1 = k + 1;
k1 = 1;
vold = zeros(n,1);
zold = zeros(n,1);
rold = zeros(n,1);
ritzold = 0;

[v,h,info,ortcost,matvectime,orttime] = arnoldno(k1,k,A,v,h,ortcost,matvectime,orttime);
matvec = matvec + k;

% compute the eigenvalues and eigenvectors of the k-by-k matrix h
[w,q1] = shftit(h,kstart,k,sigma);
ritzesone = w(k - kstart +1);

% Now we in position to update the Arnoldi decomposition
% assigning the variables
iter = 0;
Res = [];
knew = k;
ritzes = zeros(ksave,1);
ritzests = ones(ksave,1);
stopcrit = 1;
beta = 1;
betanew = 1;
residest = 1;
kkconv = 0;
kendconv = 0;
kc = zeros(2,1);
iskendconv = true;

% main loop for updating the Arnoldi sequence in place
while ((stopcrit > tol) && (iter < maxit) || (iter < 2)),
    iter = iter + 1;
    
    % computing the value of zyem
    if (kkconv == 0),
        if (iter - 1 == 0),
            [vold,zold,rold,zyem,ritzold,ritzesonesave] = ...
                zyemcomp(v,q1(:,k), ...
                ritzesone,vold,zold,rold,ritzold,iter-1,k);
        else
            [vold,zold,rold,zyem,ritzold,ritzesonesave] = ...
                zyemcomp(v,q1(:,kend), ...
                ritzesone,vold,zold,rold,ritzold,iter,kend);
        end
    else
        % setting zyem to be one if kkconv > 0 and iskendconv = false
        zyem = 1.0;
        if iskendconv,
            kendconv = kend;
            kc(1,1) = 1;
            kc(2,1) = p_max;
            iskendconv = false;
        end
    end
    
    % switching the Krylov subspace dimension
    if (iter == 1),
        [kend,kc,p_max] = switchdimen(zyem,k,k,knew, ...
            p_min,p_max,tol,kkconv,kendconv,kc,iter);
    else
        [kend,kc,p_max] = switchdimen(zyem,kend,k,knew, ...
            p_min,p_max,tol,kkconv,kendconv,kc,iter);
    end
    
    % assigning logical parameters
    isloss = true;
    isredo = false;
    
    while isloss,
        % extend the Arnoldi factorization to p additional steps
        kold = k;
        p = kend - kold;
        psave = p;
        k = knew;
        
        if isredo,
            [v,h,info,ortcost,matvectime,orttime] = arnoldno(k,kend,A,v,h,ortcost,matvectime,orttime);
        else
            [v,h,info,ortcost,matvectime,orttime] = arnoldno(k,kend,A,v,h,ortcost,matvectime,orttime);
        end
        matvec = matvec + kend - k;

        % reshaping h and v with respect to kend
        h = h(1:kend,1:kend);
        v = v(:,1:kend + 1);
        
        % setting k to the old value
        k = kold;

        % compute p shifts based on sigma to be used in implicit
        % factorization
        [w,q1] = shftit(h,kstart,kend,sigma);
        
        % sorting the ritz values with respect to sigma
        ritzes = w((kend - kstart + 1):-1:(kend - ksave + 1));
        ritzesone = ritzes(1);
             
        % checking loss of information for eigenvector
        [v,h,knew,kend,isloss,isredo,ritzesonesave] = ...
            islossevec(v,h,w,k,knew,p,ritzesone, ...
            ritzesonesave,kend,p_min,p_max,kkconv,sigma);
    end
    
    % removing display time
    cputms = cputime - t0;
    
    % updating the command window with current eigenvalue estimates
%     if display
%         ds = sprintf(['Iteration %d: a few Ritz values of the ' ...
%             '%d-by-%d matrix:'],iter,kend,kend);
%         disp(ds)
%         disp(ritzes)
%     end
    
    % restarting timing
    t0 = cputime;
    
    % computing the number of converged Ritz values
    [m1,m2] = size(q1);
%     ritz = norm(q1(m1,p + 2:m1));
    betanew = norm(v(:,kend + 1));

    ritznew = betanew*q1(m1,1:m1);
    jj = m1;
    kconv = 0;
    while (jj > 0),
        if (abs(ritznew(jj)) <= tol),
            jj = jj - 1;
            kconv = kconv + 1;
        else
            jj = -1;
        end
    end
    kkconv = kconv;

    %  At first, we use the Ritz estimates to estimate convergence.
    %  However, as the algorithm converges, the Ritz estimates become
    %  poor estimates of the actual error in each Ritz pair.  So when the
    %  Ritz estimates become too small, we actually form the matrix of
    %  errors || AV - VD || where V are the estimates for the eigenvectors
    %  and eigenvalues.  This is expensive computationally but ensures
    %  that the user gets the desired eigenpairs to the requested
    %  tolerance.
    ritzests = [w abs(ritznew)'];
    ritzests = ritzests(size(ritzests,1):-1:size(ritzests,1) - ...
        ksave + 1,2);
    Res = [Res;ritzests(1:ksave)'];
    stopcrit = max(ritzests);
    residest = norm(ritzests);
   
    if (max(ritzests) <= tol*max(norm(h),1)),
        vee = v(:,1:kend);

        errmat = A*(vee*q1(1:kend, kend:-1:kend - ksave + 1)) - ...
            vee*q1(1:kend,kend:-1:kend - ksave + 1)*diag(ritzes);
        residest = norm(errmat,1)/norm(A,1);
        matvec = matvec + ksave;

        for ii = 1:length(ritzes),
            ritzests(ii) = norm(errmat(:,ii));
        end
        stopcrit = residest;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    ve = v(:,1:kend);
    errm = A*(ve*q1(1:kend, kend:-1:kend - ksavesave + 1)) - ...
        ve*q1(1:kend,kend:-1:kend - ksavesave + 1)*diag(ritzes);
    for convi = 1:ksavesave,
        convT(iter,convi) = norm(errm(:,convi))/norm(A,1);
    end
    conv(iter,1) = matvec;
    conv(iter,2) = ortcost;
    conv(iter,3) = max(convT(iter,:));
    conv(iter,4) = matvectime;
    conv(iter,5) = orttime;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ((stopcrit > tol) || (iter < 2)),
        % Apply the p implicit shifts if convergence has not yet
        % happened.  Otherwise don't apply them and get out of the
        % loop on next loop test. We need to keep same test here as
        % in the main loop test to   avoid applying shifts and then
        % quitting, which would lead to a wrong size factorization
        % on return.

        % If some ritz values have converged then adjust k and p to
        % move the "boundary" of the filter cutoff.
        jj = 1;
        ukconv = 0;
        while (jj <= (m1 - ksave)),
            if (abs(ritznew(jj)) <= tol),
                jj = jj + 1;
                ukconv = ukconv + 1;
            else
                jj = jj + 1;
            end
        end
        kk = 0;
        if ((kconv > 0) || (ukconv > 0)),
            kk = ksavesave + kconv;
%             kk = kk + ukconv;
            p = max(ceil(psave/3),kend - kk);
            k = kend - p;
        end
        
%         % keeping the shifts
%         if (iter == 1),
%             mufilter = w(k+1:end)
%         else
%             mufilter = vertcat(mufilter,w(k+1:end))
%         end

        if (any(any(imag(v))) || any(any(imag(h)))),
            [v,h,knew] = apshft1(v,h,w,k,p);
        else
            [v,h,knew] = apshft2(v,h,real(w),imag(w),k,p);
        end

        betanew = norm(v(:,kp1));
    end
end % end of Arnoldi iteration (main loop)

if (stopcrit <= tol),
    flag = 0;
elseif (kkconv == 0)
    if (iter >= maxit),
        warning('Maximum number of iterations exceeded!')
    end

    if (nargout < 4),
        warning('None of the %d requested eigenvalues converged.',ksave)
    else
        flag = 1;
    end
elseif (kkconv < ksave),
    if (iter >= maxit),
        warning('Maximum number of iterations exceeded!')
    end
   
    if (nargout < 4),
        warning(['Only %d of the %d requested eigenvalues converged.'] ...
            , kkconv,ksave)
    else
        flag = 1;
    end
end

% compute the eigenvalues and eigenvectors of h
[w,q1] = shftit(h,1,kend,sigma);

k = ksave;
p = psave;

% transform the converged eigenvalues back to
% the original problem and return them in the diagonal
% k by k matrix d.

% set v equal to the wanted eigenvectors
v = v(:,1:kend) * q1(1:kend, kend:-1:kend - k + 1);
for i=1:k,
    v(:,i) = v(:,i)/norm(v(:,i));
end

% the eigenvalues are recovered from the set of shifts w
d = w(kend:-1:kend - k + 1);

if ~isstr(A),
    v(pmmd,:) = v;
end

if (display == 2),
%     clc
    ds = sprintf(['\n' 'Number of iteration is %d'],iter);
    disp(ds)
    
    ds = sprintf(['\n' 'The cpu timing is %d'],cputms);
    disp(ds)
    
    ds = sprintf(['\n' 'Number of matrix-vector product is %d'],matvec);
    disp(ds) 
    
    ds = sprintf(['\n' 'Number of orthogonalisation steps is %d'],ortcost);
    disp(ds)
    
    ds = sprintf(['\n' 'The maximum number of Krylov size is %d'],p_max);
    disp(ds)
end

ds = [iter matvec ortcost kkconv p_max - 1]';

% assigning the output
if (nargout <= 1),
    varargout{1} = d;
else
    varargout{1} = v;
    varargout{2} = diag(d);
    varargout{3} = ds;
    if (nargout >= 4),
        varargout{4} = conv;
%         varargout{4} = flag;
    end        
end

cputms1 = cputime - t0;
cputms = cputms + cputms1;
Res = Res/norm(A,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           arnoldno subroutine                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v,h,defloc,ortcost,matvectime,orttime] =  arnoldno(k1,k2,A,v,h,ortcost,matvectime,orttime)
% arnoldno computes or extends an Arnoldi factorization.
%
% Dan Sorensen and Richard J. Radke, 11/95.
defloc = 0;
matvectime = 0;
orttime = 0;

for j = k1 + 1:k2,
    jm1 = j - 1;
    jp1 = j + 1;
    beta = norm(v(:,j));
    h(j,jm1) = beta;
    if (beta <= 10*eps*norm(h(1:jm1,jm1))),     % If beta is "small"
        v(:,j) = rand(size(v,1),1);             % we deflate by
        s = v(:,1:jm1)'*v(:,j);                 % setting the
        v(:,j) = v(:,j) - v(:,1:jm1)*s;         % corresponding
        s = v(:,1:jm1)'*v(:,j);                 % subdiagonal of H
        v(:,j) = v(:,j) - v(:,1:jm1)*s;         % equal to 0 and
        beta = norm(v(:,j));                    % starting the
        beta = 1.0/beta;                        % basis for a new
        h(j,jm1) = 0.0;                         % invariant subspace.
        defloc = j;
    else
        beta = 1.0/beta;
    end
    
    v(:,j) = v(:,j)*beta;

    % compute w = Av and store w in j+1 -st col of V
    v2 = v(:,j);
    t = cputime;
    v(:,jp1) = A*v2;
    t = cputime - t;
    matvectime = matvectime + t;
    
    % compute the next (j-th) column of H
    t = cputime;
    h(1:j,j) = v(:,1:j)'*v(:,jp1);

    % compute the residual in the j+1 -st col of V
    v(:,jp1) = v(:,jp1) - v(:,1:j)*h(1:j,j);
    t = cputime - t;
    
    orttime = orttime + t;
    ortcost = ortcost + j;

    % perform one step of iterative refinement to correct the orthogonality
    t = cputime;
    s = v(:,1:j)'*v(:,jp1);
    v(:,jp1) = v(:,jp1) - v(:,1:j)*s;
    h(1:j,j) = h(1:j,j) + s;
    t = cputime - t;
    
    orttime = orttime + t;
    ortcost = ortcost + j;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             shftit subroutine                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,qq] = shftit(h, kstart, kend, sigma)
% shftit calculate shifts to update an Arnoldi factorization.
%
% shftit was designed to be called by dneigs.
%
% syntax: [w,q] = shftit(h, kstart, kend, sigma) where
%
% h is an upper Hessenberg matrix (from the Arnoldi factorization
%     Av = vh + fe_k'),
%
% kstart points to the start of the active block in h,
%
% kend points to the end of the active block in h,
%
% sigma is one of 'LR','SR','LM','SM'.
% 
% shftit calculates [q,w] = eig(The active block of h), where
% the eigenvalues and eigenvectors are reordered according to sigma,
% with the eigenvalues to use as shifts put first.
%
% Dan Sorensen and Richard J. Radke, 11/95.

[q,ww] = eig(h(kstart:kend,kstart:kend));
w = diag(ww);

k = kend - kstart + 1;

% select filter mechanism by activating appropriate choice below
switch (sigma),
    case {'SM'}
        % sort for smallest absolute value (shifts are largest abs val)
        [s,ir] = sort(-abs(w));
    case {'LR'}
        % sort for largest real part (shifts are smallest real part)
        [s,ir] = sort(real(w));
    case {'SR'}
        % sort for smallest real part (shifts are largest real part)
        [s,ir] = sort(-real(w));
    otherwise
        % sort for largest absolute value (shifts are smallest abs val)
        [s,ir] = sort(abs(w));
end
        
%  apply sort to w and q
w = w(ir);
qq = q(:,ir);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             apshft1 subroutine                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v,h,k] =  apshft1(v,h,w,k,p)
% apshft1 applies shifts to update an Arnoldi factorization
%
% apshft1 was designed to be called by dneigs when v or h is complex.
%
% syntax: [v,h,k] = apshft1(v,h,w,k,p) implicitly applies the
%         p real shifts held in w to update the existing Arnoldi
%         factorization Av - vh = re'_{k+p}.
%                               
% the routine results in
%
%     A(vq) - (vq)(q'hq = re'_{k+p} *q
%                            
% where the orthogonal matrix q is the product of the Givens
% rotations resulting from p bulge chase sweeps.
%
% the updated residual is placed in v(:,k+1) and the updated
% Arnoldi factorization v <- vq, h <- q'hq is returned.
%
% Dan Sorensen and Richard J. Radke, 11/95.

k1 = 1;
kend = k + p;
kp1 = k + 1;
q = zeros(1,kend);
q(kend) = 1.0;

ix = find(diag(h,-1) == 0);     % Find the column indices of 0
ix = [0 ; ix ; kend - 1];       % subdiagonals in h.
nx = size(ix,1);

for jj = 1:p,               % Loop over shifts
    for ii = 1:nx - 1,         % Loop over blocks in h
        k1 = ix(ii) + 1; k2 = ix(ii + 1);
        c = h(k1,k1) - w(jj);
        s = h(k1 + 1,k1);
        
        [G,R] = qr([c;s]);
        for i = k1:k2,        % Loop over rows in the block
            if (i > k1),
                [G,R] = qr(h(i:i + 1,i - 1));
                h(i:i + 1,i - 1) = R(:,1);
            end

            % apply rotation from left to rows of h
            h(i:i + 1,i:kend) = G'*h(i:i + 1,i:kend);

            % apply rotation from right to columns of h
            ip2 = i + 2;
            if (ip2 > kend),
                ip2 = kend;
            end
            h(1:ip2,i:i + 1) = h(1:ip2,i:i + 1)*G;

            % apply rotation from right to columns of v
            v(:,i:i + 1) = v(:,i:i + 1)*G;

            % accumulate e'_{k+p} *q so residual can be updated
            q(i:i + 1) = q(i:i + 1)*G;
        end
    end
end

% update the residual and store in the k+1 -st column of v
v(:,kend + 1) = v(:,kend + 1)*q(k);
v(:,kp1) = v(:,kend + 1) + v(:,kp1)*h(kp1,k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             apshft2 subroutine                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v,h,k] =  apshft2(v,h,wr,wi,k,p)
% apshft2 applies shifts to update an Arnoldi factorization
%
% apshft2 was designed to be called by dneigs when V and H are real.
%
% syntax: [v,h,k] = apshft2(v,h,wr,wi,k,p) implicitly applies
%     the p complex shifts given by w = wr + i*wi, to update the
%     existing Arnoldi factorization Av - vh = re'_{k+p}.
%
% the routine results in
%
%     A(vq) - (vq)(q'hq) = re'_{k+p} *q
%
% where the orthogonal matrix q is the product of the Givens
% rotations resulting from p bulge chase sweeps.
%
% The updated residual is placed in v(:,k+1) and the updated
% Arnoldi factorization v <- vq, h <- q'hq is returned.
%
% Dan Sorensen and Richard J. Radke, 11/95.

k1 = 1;
kend = k + p;
kp1 = k + 1;
q = zeros(kend,1);
q(kend) = 1.0;
mark = 0;
num = 0;

ix1 = find(diag(h,-1) == 0);    % Find the column indices of 0
ix = [0 ; ix1 ; kend - 1];         % subdiagonals in h.
jx = [0 ; ix1 ; kend - 2];
nx = size(ix,1);

for jj = 1:p,
    % compute and apply a bulge chase sweep initiated by the
    % implicit shift held in w(jj)
    if (abs(wi(jj)) == 0.0),
        % apply a real shift using 2 by 2 Givens rotations
        for ii = 1:nx - 1,      % Loop over blocks in h
            k1 = ix(ii) + 1;
            k2 = ix(ii + 1);
            
            c = h(k1,k1) - wr(jj);
            s = h(k1 + 1,k1) ;
            t = norm([c s]);
            
            if (t == 0.0),
                c = 1.0;
            else
                c = c/t;
                s = s/t;
            end
            
            for i = k1:k2,
                if (i > k1),
                    t = norm(h(i:i + 1,i - 1));
                    if (t == 0.0),
                        c = 1.0;
                        s = 0.0;
                    else
                        c = h(i,i - 1)/t;
                        s = h(i + 1,i - 1)/t;
                        h(i,i - 1) = t;
                        h(i + 1,i - 1) = 0.0;
                    end
                end

                % apply rotation from left to rows of h
                G = [c s ; -s c];
                h(i:i + 1,i:kend) = G* h(i:i + 1,i:kend);

                % apply rotation from right to columns of h
                ip2 = i + 2;
                
                if (ip2 > kend),
                    ip2 = kend;
                end
                h(1:ip2,i:i + 1) =  h(1:ip2,i:i + 1)*G';

                % apply rotation from right to columns of v
                v(:,i:i + 1) =  v(:,i:i + 1)*G';

                %  accumulate e'_{k+p} *q so residual can be updated            
                q(i:i + 1) =  G*q(i:i + 1);
            end
        end
        num = num + 1;
    else
        % apply a double complex shift using 3 by 3 Householder
        % transformations
        if ((jj == p) || (mark == 1)),
            mark = 0;       % skip application of conjugate shift
        else
            num = num + 2;    % mark that a complex conjugate
            mark = 1;         % pair has been applied

            for ii = 1:nx - 1,   % Loop over blocks in H
                k1 = jx(ii) + 1;
                k2 = k1 + 1;
                k3 = jx(ii + 1);
                c = h(k1,k1)*h(k1,k1) + h(k1,k2)*h(k2,k1) ...
                    - 2.0*wr(jj)*h(k1,k1);
                c = c + wr(jj)^2 + wi(jj)^2;
                s = h(k2,k1)*(h(k1,k1) + h(k2,k2) - 2.0*wr(jj));
                g = h(k2+1,k2)*h(k2,k1);
                t = norm([c s g]);
                sig = -sign(c);
                c = c -t*sig;
                for i = k1:k3,
                    if (i > k1),
                        t = norm(h(i:i + 2,i - 1));
                        sig = -sign(h(i,i - 1));
                        c = h(i,i - 1) - t*sig;
                        s = h(i + 1,i - 1);
                        g = h(i + 2,i - 1);
                        h(i,i - 1) = t;
                        h(i + 1,i - 1) = 0.0;
                        h(i + 2,i - 1) = 0.0;
                    end
                    t = norm([c s g]);
                    if (t ~= 0.0),
                        c = c/t;
                        s = s/t;
                        g = g/t;
                    end
                    z = [c s g]';

                    % apply transformation from left to rows of h
                    t =  sig*2.0*(z'*h(i:i + 2,i:kend));
                    h(i:i + 2,i:kend) = sig*h(i:i + 2,i:kend) - z*t;

                    % apply transformation from right to columns of h
                    ip3 = i + 3;
                    if (ip3 > kend),
                        ip3 = kend;
                    end

                    t =  sig*2.0*h(1:ip3,i:i + 2)*z;
                    h(1:ip3,i:i + 2) = sig*h(1:ip3,i:i + 2) - t*z';

                    % apply transformation from right to columns of v
                    t =  sig*2.0*v(:,i:i + 2)*z;
                    v(:,i:i + 2) = sig*v(:,i:i + 2) - t*z';

                    % accumulate e'_{k+p} *q so residual can be updated
                    t =  sig*2.0*z'*q(i:i + 2);
                    q(i:i + 2) = sig*q(i:i + 2) - z*t;
                end
            end

            % clean up step with Givens rotation
            i = kend - 1;
            if (i > k1),
                t = norm([h(i,i - 1) h(i + 1,i - 1)]);
                if (t ~= 0.0),
                    c = h(i,i - 1)/t;
                    s = h(i + 1,i - 1)/t;
                else
                    c = 1.0;
                    s = 0.0;
                end
                h(i,i - 1) = t;
                h(i + 1,i - 1) = 0.0;
            end

            % apply rotation from left to rows of h
            G = [c s ; -s c];
            h(i:i + 1,i:kend) = G* h(i:i + 1,i:kend);

            % apply rotation from right to columns of h
            ip2 = i + 2;
            if (ip2 > kend),
                ip2 = kend;
            end
            h(1:ip2,i:i + 1) =  h(1:ip2,i:i + 1)*G';

            % apply rotation from right to columns of v
            v(:,i:i + 1) =  v(:,i:i + 1)*G';

            % accumulate e'_{k+p} *q so residual can be updated            
            q(i:i + 1) =  G*q(i:i + 1);
        end
    end
end

% update residual and store in the k+1 -st column of v
k = kend - num;
v(:,kend + 1) = v(:,kend + 1)*q(k);

if (k < size(h,1)),
    v(:,k + 1) = v(:,kend + 1) + v(:,k + 1)*h(k + 1,k);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            zyemcomp subroutine                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vold,zold,rold,zyem,ritzold,ritzesonesave] = zyemcomp(v,q1, ...
            ritzesone,vold,zold,rold,ritzold,itr,kend)
% zyemcomp was designed to be called by dneigs.
%
% zyemcomp computes the value of zyem that will be used to switch the
% Krylov subspace dimension dynamically.
%
% Dookhitram Kumar and Bhuruth Muddun, 12/2007.

if (itr == 0),
    % computing the eigenvector corresponding to the first eigenvalue
    % respective to sigma
    x = v(:,1:kend)*q1;
    beta = 1.0/norm(x);
    x = x*beta;
    
    % computing the residual vector with respect to x
    r = q1(kend)*v(:,kend+1);
    
    % adding the residual vector with a multiple of x
    z = ritzesone*x + r;
    
    % set zyem to be one to start with a small Krylov subspace dimension
    zyem = 1.0;
    
    % replacing the value of old variables with new value
    vold = x;
    zold = z;
    rold = r;
    ritzold = ritzesone;
    
    % saving the value of ritzesone
    ritzesonesave = ritzesone;
else
    % computing the eigenvector corresponding to the first eigenvalue
    % respective to sigma
    x = v(:,1:kend)*q1;
    beta = 1.0/norm(x);
    x = x*beta;
    
    % computing the residual vector with respect to x
    r = q1(kend)*v(:,kend+1);
    
    % adding the residual vector with a multiple of x
    z = ritzesone*x + r;
    
    % compting the value of eta
    eta = z - zold - ritzold*(x - vold) - (ritzesone - ritzold)*x;
    
    % evaluating zyem using rold and eta
    zyem = rold'*eta;norm(rold);
    beta = 1.0/(norm(rold)*norm(eta));
    zyem = zyem*beta;
    
    % recomputing r by using rold and eta
    r = rold + eta;
    
    % replacing the value of old variables with new value
    vold = x;
    zold = z;
    rold = r;
    ritzold = ritzesone;
    
    % saving the value of ritzesone
    ritzesonesave = ritzesone;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           switchdimen subroutine                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [kend,kc,p_max] = switchdimen(zyem,kend,k, ...
    knew,p_min,p_max,tol,kkconv,kendconv,kc,iter)
% switchdimen was designed to be called by dneigs.
%
% switchdimen switches the Krylov subspace dimension with respect to the
% value of zyem.
%
% Dookhitram Kumar and Bhuruth Muddun, 12/2007. 

if ((k - knew) <= 0),
    kswitch = knew;
else
    kswitch = k;
end

if kc(1,1),
    kend = kc(2,1);
    kc(1,1) = 0;
end

if (abs(zyem) < (1.0 - tol)),
    if (k > 7),
        if (mod(iter,6) == 0),
            kend = kend + 0;
        else
            kend = kend + 1;
        end
        if (kend >= p_max),
            kend = p_min;
        end
    else
        if (mod(iter,6) == 0),
            kend = kend - 0;
        else
            kend = kend - 1;
        end
        if (kend < p_min),
            kend = p_max - 1;
        end
    end
elseif ((abs(zyem) > (1.0 - tol)) && (abs(zyem) <= 1)),
    if (kkconv == 0),
        if (k > 7),
            kend = p_min;
        else
            kend = p_max - 1;
        end
    else
        if (((kend - 1) > p_min) && ((kend - kswitch - 1) >= 2)),
            if (k > 7),
                kend = kend - 1;
            else
                if (mod(iter,6) == 0),
                    kend = kend - 1;
                else
                    kend = kend - 3;
                end
            end
        else
            if ((p_max - kswitch - 1) >= 2),
                kend = p_max - 1;
            elseif ((p_max - kswitch - 1) == 1),
                p_max = p_max + 1;
                kend = p_max - 1;
            elseif ((p_max - kswitch - 1) == 0),
                p_max = p_max + 2;
                kend = p_max - 1;
            end
        end
    end
end

% if (abs(zyem) < (1.0 - tol)),
%     kend = kend - 1;
%     if (kend <= p_min),
%         kend = p_max - 1;
%     end
% elseif ((abs(zyem) > (1.0 - tol)) && (abs(zyem) <= 1)),
%     if (kkconv == 0),
%         kend = p_max - 1;
%     else
%         kend = kend - 1;
%         if ((kend <= kendconv) && (kend <= kswitch)),  
%             kend = p_max - 1;
%         end
%     end
% end
% if (abs(zyem) < (1.0 - tol)),
%     kend = kend - 1;
%     if (kend <= p_min),
%         kend = p_max - 1;
%     end
% elseif ((abs(zyem) > (1.0 - tol)) && (abs(zyem) <= 1)),
%     if (kkconv == 0),
%         kend = p_max - 1;
%     else
%         kend = kend + 1;
%         if ((kend >= p_max) && ((kswitch + 1) < (p_max - 1))),
%             if ((kendconv > kswitch) && (kendconv ~= (p_max - 1))),
%                 kend = kendconv;
%             else
%                 kend = kswitch + 1;
%             end
%         else
%             kend = p_max - 1;
%         end
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          islossevec  subroutine                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v,h,knew,kend,isloss,isredo,ritzesonesave] = ...
    islossevec(v,h,w,k,knew,p,ritzesone,ritzesonesave, ...
    kend,p_min,p_max,kkconv,sigma);
% islossevec was designed to be called by dneigs.
%
% islossevec check if there is loss of information for eigenvector
% corresponding to sigma
%
% Dookhitram Kumar and Bhuruth Muddun, 12/2007.

if (kkconv == 0),
    x = kend + fix(k/3);
    switch (sigma),
        case {'LR'}
            if ((real(ritzesonesave) >= real(ritzesone)) && (x < p_max)),
                if (any(any(imag(v))) || any(any(imag(h)))),
                    [v,h,knew] = apshft1(v,h,w,k,p);
                else
                    [v,h,knew] = apshft2(v,h,real(w),imag(w),k,p);
                end
                kend = x;
                isloss = true;
                isredo = true;
            else
                ritzesonesave = ritzesone;
                isloss = false;
                isredo = false;
            end
        case {'SR'}
            if ((real(ritzesonesave) <= real(ritzesone)) && (x < p_max)),
                if (any(any(imag(v))) || any(any(imag(h)))),
                    [v,h,knew] = apshft1(v,h,w,k,p);
                else
                    [v,h,knew] = apshft2(v,h,real(w),imag(w),k,p);
                end
                kend = x;
                isloss = true;
                isredo = true;
            else
                ritzesonesave = ritzesone;
                isloss = false;
                isredo = false;
            end
        case {'SM'}
            if ((abs(ritzesonesave) <= abs(ritzesone)) && (x < p_max)),
                if (any(any(imag(v))) || any(any(imag(h)))),
                    [v,h,knew] = apshft1(v,h,w,k,p);
                else
                    [v,h,knew] = apshft2(v,h,real(w),imag(w),k,p);
                end
                kend = x;
                isloss = true;
                isredo = true;
            else
                ritzesonesave = ritzesone;
                isloss = false;
                isredo = false;
            end
        otherwise
            if ((abs(ritzesonesave) >= abs(ritzesone)) && (x < p_max)),
                if (any(any(imag(v))) || any(any(imag(h)))),
                    [v,h,knew] = apshft1(v,h,w,k,p);
                else
                    [v,h,knew] = apshft2(v,h,real(w),imag(w),k,p);
                end
                kend = x;
                isloss = true;
                isredo = true;
            else
                ritzesonesave = ritzesone;
                isloss = false;
                isredo = false;
            end
    end
else
    ritzesonesave = ritzesone;
    isloss = false;
    isredo = false;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%