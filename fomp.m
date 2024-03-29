% Fast omp test.
% Really just testing Householder form of dense QR.
% Santiago Akle, ICME, Stanford University, tiagoakle@gmail.com
% 26 Nov 2012: Polished original version.

clear all
% Set the experiment parameters
m   = 16;
n   = 32;
k   = 2;
tol = 1e-10;

% Select a support
rand('seed',0);
randn('seed',0);
T = randperm(n);
T = T(1:k);

% Select a random measurement set
M = randperm(n);
M = M(1:m);
M = sort(M);

% Generate a sparse solution
x_s    = zeros(n,1);
x_s(T) = randn(k,1);

% Generate the measurement matrix
A      = zeros(m,n);
e      = zeros(1,n);
for i=1:m
   e(M(i)) = 1;
   A(i,:)  = dct(e);
   e(M(i)) = 0;
end
% Form the measurement vector
b      = A*x_s;

% Set up the omp variables
S      = zeros(m,k);
S_i    = zeros(k,1); % Indices vector
Q_v    = zeros(m,m); % Matrix of householder vectors
b_h    = b;          % Rotated b vector
r      = b;          % Present residual
x      = [];         % Present solution
bnorm  = norm(b);    % Residual norm
rnorm  = bnorm;

% Iteration variables and other stuff
itn    = 1;          % iteration count
opts.UT= true;       % linsolve options
e       = zeros(n,1);

% Main OMP iteration
while itn<m && rnorm/bnorm > tol
    z          = A'*r;        % We want to substitute this for sfft
    [mp,p]     = max(abs(z)); % This is an order n operation if we do it this way but sfft should improve
    S_i(itn)   = p;           % Save the index of the maximul col 

    % Extract the column of A
    a          = A(:,p); 
    S(:,itn)   = a;
    if itn>1
        a      = product_q(Q_v,2,a); % Apply the Householder rotation in product form to this new column 
                                     % (requires itn inner products of size (m-itn+1) in the last iteration O(mk) )
    end
    v              = house(a(itn:m));      % Calculate the new reflection vector. This operation is O(m-itn+1) (calculate a norm)
    t              = v'*a(itn:m);          % Apply the rotation to the rest of a to form the new column of R
    a(itn)         = a(itn) - 2*v(1)*t;    
    Q_v(itn:m,itn) = v;                    % Save the new vector
    R(1:itn,itn)   = a(1:itn);             % Save the new column of R
    t              = v'*b_h(itn:m);
    b_h(itn:m)     = b_h(itn:m) - (2*t)*v; % Apply the new rotation to b_h O(m-itn)

    x     = linsolve(R,b_h(1:itn),opts);   % Solve the least squares problem O(itn^2)
    
    r     = b - S(:,1:itn)*x;              % Calculate the new residual, all indices changed so it is not that cheap O(m*itn);
    rnorm = norm(r);                       % Norm for the termination criteria O(m)
    itn   = itn + 1;                                         
end
itn = itn - 1;
