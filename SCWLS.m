function [funML, x_est, times] = SCWLS(array, realization, Q_mat)
% new Constarined WLS algorithm for TDOA-based positioning (proposed)
% The new twist here is that G'G is nosigular now when unknown-position
% sensor is located in the certral of the anchors
% in general cwls method, A = [G g]; theta = [eta R1]
% model: \A \eta = \b - \g R_1+ \m ; subject to \eta^T \eta = R_1^2 
% where : \eta = [x-x_1 y-y_1]^T                  
% -----------------------------------------------------------------------
% x = cwls_tdoa_proposed(X,r,sigma2)
% X = Anchors position
% x = 2D position estimate
% r_i1 = TDOA measurement vector, r_i1 = r_i - r_1, i=2,3,...,M

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use the standard root finding technique in Gary's CWLS algorithm to find
% lambda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% 2023 update %%%%%
SensorPositions = realization.sensors;
r_i1 = realization.d(2:end)';
N = array.N;
n = array.n;
d = realization.d;
array_no_ref = SensorPositions(:, 2:N); ref = SensorPositions(:, 1);
ML = array.ML;
if contains(Q_mat, 'Id')
    Q = eye(N - 1);
elseif contains(Q_mat, 'sig')
    sigma = realization.sigma;
    Q = (sigma^2)*ones(N - 1) + (sigma^2)*eye(N - 1);
end
SensorPositions = SensorPositions';
%%%%%%%%%%%%%%%%%%%%%%%

% [N dim] = size(X);
X_i1 = SensorPositions(2:N,:)-ones(N-1,1)*SensorPositions(1,:);
G = 2*X_i1;
g = 2*r_i1;
b = sum(X_i1.^2,2)-r_i1.^2;

% sigma2_i1= sigma2(1)+sigma2(2:end);
D = eye(N-1);
% Q = diag(sigma2_i1)+(ones(N-1)-eye(N-1))*sigma2(1);
C = D*Q*D;

Lw = 3;

A = [G g];
S = diag([1 1 -1]);

tic
for index_w = 1:Lw
    
    W = inv(C);
    [U,L] = eig(A'*W*A*S);
    pp = U'*S*A'*W*b;
    qq = (U)\(A'*W*b);
    p1 = pp(1);
    p2 = pp(2);
    p3 = pp(3);
    q1 = qq(1);
    q2 = qq(2);
    q3 = qq(3);
    r1 = L(1,1);
    r2 = L(2,2);
    r3 = L(3,3);
    coeff1=p1*q1+p2*q2+p3*q3;
    coeff2=2*p2*q2*r1+2*p3*q3*r1+2*p1*q1*r2+2*p3*q3*r2+2*p1*q1*r3+2*p2*q2*r3;
    coeff3=p2*q2*r1^2+p3*q3*r1^2+4*p3*q3*r1*r2+p1*q1*r2^2+p3*q3*r2^2+4*p2*q2*r1*r3+4*p1*q1*r2*r3+p1*q1*r3^2+p2*q2*r3^2;
    coeff4=2*p3*q3*r1^2*r2+2*p3*q3*r1*r2^2+2*p2*q2*r1^2*r3+2*p1*q1*r2^2*r3+2*p2*q2*r1*r3^2+2*p1*q1*r2*r3^2;
    coeff5=p3*q3*r1^2*r2^2+p2*q2*r1^2*r3^2+p1*q1*r2^2*r3^2;

    polynomial = [coeff1 coeff2 coeff3 coeff4 coeff5];

    root = roots(polynomial);

    %***********find the smallest root(lambda)********************** 
    clear storage_root IndexJn;
    clear Jn;
    counter = 0;
    for i = 1:length(root)
        if isreal(root(i)) 
           counter = counter + 1; 
           storage_root(counter) = root(i);
        end
    end
    Jn = zeros(counter,1);
    R1_est = zeros(counter,1);
    eta_est = zeros(n,counter);
    
    num = zeros(N,1);

    for i=1:counter
        [R1_est(i), num(i)] = R_root_selection_tdoa(G,g,b,W,storage_root(i),num(i));
        eta_est(:,i) = (G'*W*G+storage_root(i)*eye(n))\(G'*W*b-G'*W*g*R1_est(i));
        Jn(i) = (G*eta_est(:,i)-b+g*R1_est(i))'*W*(G*eta_est(:,i)-b+g*R1_est(i));

    end
    [~, IndexJn]=sort(Jn);
    % lambda=storage_root(IndexJn(1));

    eta = eta_est(:,IndexJn(1));
    x_est = eta'+SensorPositions(1,:);
    d_est = sqrt(sum((SensorPositions - ones(N,1)*x_est).^2,2));
    D = 2*diag(d_est(2:end));
    C = D*Q*D;

end
times = toc;

%cd = cond(G'*W*G+lambda*eye(dim));

%%%%% 2023 update %%%%%
funML = ML(x_est',d,array_no_ref,ref);
%%%%%%%%%%%%%%%%%%%%%%%
end



function [R1 num] = R_root_selection_tdoa(G,g,b,W,lambda,num)
% Root Selection Procedure for R1
% R1 = sqrt((x-x_1)^2+ (y-y1)^2)
% To solve the equation: a1*R1^2+b1*R1+c1=0
% ---------------------------------------------
% R1 = R_root_selection_tdoa(G,g,b,W,lambda)
% Cost function:J = (G*eta-b+g*R1)^T*W*(G*eta-b+g*R1)
% source position: [x y]^T
% ith receiver position: [x_i y_i]^T
% eta = [x-x1; y-y1];
% d_i = sqrt((x-x_i)^2+(y-y_i)^2),distance between source and ith receiver
% r_{i,1} = r_i -r_1; and r_i = d_i+n_i
% G = [x_2-x_1 y_2 -y_1; .... ;x_M-x_1 y_M -y_1]
% g = [r_{2,1} .... r_{M,1}]^T;
% b = 0.5*[(x_2-x_1)^2+(y_2-y_1)^2-r_{2,1}^2 .....
% (x_M-x_1)^2+(y_M-y_1)^2-r_{M,1}^2]^T
% ----------------------------------------------

dim = size(G,2);
t = inv(G'*W*G+lambda*eye(dim)); 
T = t*t;
t1 = G'*W*g;
t2 = G'*W*b;
a1 = t1'*T*t1-1;
b1 = -2*t1'*T*t2;
c1 = t2'*T*t2;

delta = b1^2-4*a1*c1;
num = 0;
if delta>0
    delta_r = sqrt(delta);
    R1_est1(1)=(-b1+delta_r)/(2*a1);
    R1_est1(2)=(-b1-delta_r)/(2*a1);
    if R1_est1(1)>=0 && R1_est1(2)<0
        R1 = R1_est1(1);
    else if R1_est1(2)>=0 && R1_est1(1)<0
            R1 = R1_est1(2);
        else if R1_est1(2)>=0 && R1_est1(1)>=0
                eta(:,1) = t*(t2 - t1*R1_est1(1));
                K1 = G*eta(:,1)-b+g*R1_est1(1);
                JJ(1) = K1'*W*K1;
                eta(:,2) = t*(t2 - t1*R1_est1(2));
                K2 = G*eta(:,2)-b+g*R1_est1(2);
                JJ(2) = K2'*W*K2;

                if JJ(1)<JJ(2)
                    R1 = R1_est1(1);
                else
                    R1 = R1_est1(2);
                end

            else
                R1 = 0;
                num = 1;
            end
        end
    end

else
    num = 1;
    if -b1/(2*a1)>=0
        R1=-b1/(2*a1);
    else
        R1=0;
    end
end
end