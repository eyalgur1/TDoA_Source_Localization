function [funML, u, times] = SDP(array, realization, eta)

S = realization.sensors;
t = realization.d(2:end)';
N = array.N;
n = array.n;
dist = realization.d;
array_no_ref = S(:, 2:N); ref = S(:, 1);
ML = array.ML;

A = [-ones(N-1,1) eye(N-1)];

tic
cvx_begin sdp
variables u(n) d(N) Y
variable D(N,N) symmetric 
minimize square_pos(norm(A*d)) - 2*(t')*A*d + eta*square_pos(norm(D,'fro'))
subject to
for i=1:N
norm(u-S(:,i)) <= d(i)
D(i,i) == Y - 2*u'*S(:,i) + norm(S(:,i))^2
end
for i=1:N-1
    for j=i+1:N
        D(i,j) >= norm(Y - u'*S(:,i) - u'*S(:,j) + S(:,i)'*S(:,j),Inf)
    end
end
Y >= square_pos(norm(u))
[1 d';d D] == semidefinite(N+1);
cvx_end
times = toc;

funML = ML(u,dist,array_no_ref,ref);