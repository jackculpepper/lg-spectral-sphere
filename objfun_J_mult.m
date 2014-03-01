function [f,g] = objfun_J(x0,t,I0,I1,gamma);

[L,batch_size] = size(I0);
M = size(t,1);

if batch_size > 1
    fprintf('batch_size > 1 not supported yet..\n');
    keyboard
end

J = reshape(x0, L, L, M);

A = zeros(L);
for i = 1:M
    A = A + J(:,:,i) * t(i);
end

tic

[U,D] = eig(A);
V = inv(U).';

D = diag(D);

time_eigs = toc;


tic
ExpA = U * diag(exp(D)) * V.';

%f = trace(I1'*ExpA*I0);
f = -2 * trace(I1'*ExpA*I0) + trace(I0'*ExpA'*ExpA*I0) + 0.5*gamma*sum(J(:).^2);
time_obj = toc;


tic

%% populate F
F = zeros(L,L);

ExpD = exp(D);
for i = 1:L
    for j = 1:L
        if D(i) == D(j)
            F(i,j) = ExpD(i);
        else
            F(i,j) = (ExpD(j) - ExpD(i)) / (D(j) - D(i));
        end
    end
end

time_popF = toc;


tic

B = I1*I0';
Q = U.' * B * V;
R = F .* Q;
S = V * R * U.';

dJ1 = -2*S(:) * t';


B = (I0*I0'*ExpA')';
Q = U.' * B * V;
R = F .* Q;
S = V * R * U.';

dJ2 = 2*S(:) * t';

dJ = dJ1 + dJ2 + gamma*reshape(J,L^2,M);

time_mult = toc;

%fprintf('obj %.4f eigs %.4f popF %.4f mult %.4f\n', time_obj, ...
%    time_eigs, time_popF, time_mult);

g = dJ(:);

%keyboard ; barf

