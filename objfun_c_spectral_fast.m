function [f,g] = objfun_c(x0,J,I0,I1,lambda);

[L,batch_size] = size(I0);
Lsz = sqrt(L);

c = x0;

M = size(c, 1);

A = zeros(L);
for i = 1:M
    A = A + J(:,:,i) * c(i);
end

[U,D] = eig(A);
V = inv(U).';
D = diag(D);

ExpA = real(expm(A));



E = I1 - ExpA*I0;

%if (rand < 0.25 && ceil(Lsz) == Lsz)
%
%    %cvals = '';
%    %for i = 1:M
%    %    cvals = sprintf('%s c%d %.4f', cvals, i, c(i));
%    %end
%
%    figure(5);
%    subplot(2,2,3); imagesc(reshape(ExpA*I0, Lsz, Lsz)); colormap(gray);
%        axis image; colorbar;
%    subplot(2,2,4); imagesc(reshape(E, Lsz, Lsz)); colormap(gray);
%        axis image; title(sprintf('%.4f', sum(E(:).^2))); colorbar;
%
%    figure(6);
%    bar(c);
%
%    drawnow;
%end





%f = trace(I1'*ExpA*I0);
f1 = -2 * trace(I1'*ExpA*I0);
f2 = trace(I0'*ExpA'*ExpA*I0);
%f3 = 0.5*lambda*sum(c(:).^2);
f3 = 0.5*lambda*sum(abs(c(:)));

fprintf('\rf1 %.4f f2 %.4f f3 %.4f ', f1, f2, f3);

f = f1+f2+f3;



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


J_2d = reshape(J, L^2, M);

B = I1*I0';
Q = U.' * B * V;
R = F .* Q;
S = V * R * U.';

dc1 = -2*S(:)' * J_2d;



B = (I0*I0'*ExpA')';
Q = U.' * B * V;
R = F .* Q;
S = V * R * U.';

dc2 = 2*S(:)' * J_2d;

dc = dc1 + dc2 + lambda*sign(c)';

%if (rand < 0.25 && ceil(Lsz) == Lsz)
%    figure(7);
%    bar(dc);
%    drawnow;
%end


g = dc(:);

