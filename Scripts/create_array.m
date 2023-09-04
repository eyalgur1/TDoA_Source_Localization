function array = create_array(N, n, positions, s_real)

I = eye(n); e1 = I(:,1); E1 = [1 -ones(1, N-1)];

array.N = N;
array.n = n;
array.I = I;
array.e1 = e1;
array.E1 = E1;
array.s_real = s_real;
array.positions = positions;


% function handles
phi = @(s,u,p,d,array_no_ref,ref)(N-1)*norm(s)^2 - p'*s + ...
    sum((sum(u(:,2:N).*(s-array_no_ref)))...
    .*(u(:,1)'*(s-ref) + d(2:N)));  % the smooth function phi(s,u)

psi = @(s,d,ref)sum(d(2:N).*norm(s-ref));  % the non-smooth function psi(s)

Psi = @(s,u,p,d,array_no_ref,ref)phi(s,u,p,d,array_no_ref,ref) + ...
    psi(s,d,ref);  % the reformulation function Psi(s,u)

ML = @(s,d,array_no_ref,ref)sum((sqrt(sum((s-array_no_ref).^2))'...
    -d(2:N)'-norm(s-ref)).^2);  % ML formulation

array.phi = phi;
array.psi = psi;
array.Psi = Psi;
array.ML = ML;