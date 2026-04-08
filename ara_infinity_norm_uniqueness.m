clear
clc;

rng(1);

N = 100;
n = 5;
m = 2;
    
X = randn(N,n); % Data matrix
V = randn(n,m); % Matrix of axis vectors

X = zscore(X); % Standardize data


%% ARA unconstrained with the l-infinity norm

% First ARA formulation 
cvx_begin

    cvx_quiet(true);    
    variable P(N,m);
    variable t(N,1);

    minimize( max(t) );

    P*V' - X <= t*ones(1,n)
    -t*ones(1,n) <= P*V' - X

cvx_end 
P1 = P;

z1_vec = t;
z1 = max(z1_vec);  % Optimal value


% Second ARA formulation 
cvx_begin

    cvx_quiet(true);    
    variable P(N,m);
    variable t(N,1);

    minimize( sum(t) );

    P*V' - X <= t*ones(1,n)
    -t*ones(1,n) <= P*V' - X

cvx_end 
P2 = P;

z2_vec = t;
z2 = max(z2_vec);  % Optimal value


% The solutions P1 and P2 to both problems are different, but yield the
% same objective value z1==z2

% Thus, we solve the quadratic problem:
cvx_begin
    cvx_quiet(true);    
    variable P(N,m);

    minimize( norm(P*V' - X,'fro') );

    -z1*ones(N,n) <= P*V' - X
    P*V' - X <= z1*ones(N,n)

cvx_end 
P_extended = P;

% The original objective value does not change (z==z1==z2)
z = max(max(abs(P_extended*V'-X)));  

% But the associated l-2 norm is improved (val_extended<val1, and val_extended<val_2)
val_1 = norm(P1*V'-X,'fro');
val_2 = norm(P2*V'-X,'fro');
val_extended = norm(P_extended*V'-X,'fro'); 

 


%% ARA exact with the l-infinity norm

k = 1; % Selected variable

v = V(k,:)';  % k-th axis vector
x = X(:,k);   % k-th data column 

% First ARA formulation 
cvx_begin

    cvx_quiet(true);    
    variable P(N,m);
    variable t(N,1);

    minimize( max(t) );

    P*V' - X <= t*ones(1,n)
    -t*ones(1,n) <= P*V' - X

    P*v == x

cvx_end 
P1 = P;

z1_vec = t;
z1 = max(z1_vec);


% Second ARA formulation 
cvx_begin

    cvx_quiet(true);    
    variable P(N,m);
    variable t(N,1);

    minimize( sum(t) );

    P*V' - X <= t*ones(1,n)
    -t*ones(1,n) <= P*V' - X

    P*v == x

cvx_end 
P2 = P;

z2_vec = t;
z2 = max(z2_vec);


% Quadratic problem:
cvx_begin
    cvx_quiet(true);    
    variable P(N,m);

    minimize( norm(P*V' - X,'fro') );

    -z1*ones(N,n) <= P*V' - X
    P*V' - X <= z1*ones(N,n)

    P*v == x

cvx_end 
toc
P_extended = P;


% The original objective value does not change (z==z1==z2)
z = max(max(abs(P_extended*V'-X)));  

% But the associated l-2 norm is improved (val_extended<val1, and val_extended<val_2)
val_1 = norm(P1*V'-X,'fro');
val_2 = norm(P2*V'-X,'fro');
val_extended = norm(P_extended*V'-X,'fro'); 



%% ARA ordered with the l-infinity norm

[~, I] = sort(X(:,k));  % sort indices for k-th variable

% First ARA formulation 
cvx_begin

    cvx_quiet(true);    
    variable P(N,m);
    variable t;

    minimize( t );

    P*V' - X <= t*ones(N,n)
    -t*ones(N,n) <= P*V' - X

    P(I(1:N-1),:)*V(k,:)' <= P(I(2:N),:)*V(k,:)'

cvx_end 
P1 = P;

z1 = t;


% Second ARA formulation 
cvx_begin

    cvx_quiet(true);    
    variable P(N,m);
    variable t(N,1);

    minimize( max(t) );

    P*V' - X <= t*ones(1,n)
    -t*ones(1,n) <= P*V' - X

    P(I(1:N-1),:)*V(k,:)' <= P(I(2:N),:)*V(k,:)'

cvx_end 
P2 = P;

z2_vec = t;
z2 = max(z2_vec);


% Quadratic problem:
cvx_begin
    cvx_quiet(true);    
    variable P(N,m);

    minimize( norm(P*V' - X,'fro') );

    P(I(1:N-1),:)*V(k,:)' <= P(I(2:N),:)*V(k,:)'

    -z1*ones(N,n) <= P*V' - X
    P*V' - X <= z1*ones(N,n)

cvx_end 
toc
P_extended = P;


% The original objective value does not change (z==z1==z2)
z = max(max(abs(P_extended*V'-X)));  

% But the associated l-2 norm is improved (val_extended<val1, and val_extended<val_2)
val_1 = norm(P1*V'-X,'fro');
val_2 = norm(P2*V'-X,'fro');
val_extended = norm(P_extended*V'-X,'fro'); 



