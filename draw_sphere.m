load state/L=003_20090513T221016/J.mat
L = 3;
M = 9;
Msz = sqrt(M);

K(:,:,1) = J(:,:,4);
K(:,:,2) = J(:,:,6);
K(:,:,3) = J(:,:,9);
J = K;

figure(8); colormap(gray);
for i = 1:size(J,2)
    subplot(size(J,2),1,i); imagesc(J(:,:,i)); axis image; colorbar;
end

range = -10000:10:10000;
range = -200:1:200;

figure(4);
clf;

colors = { 'b', 'r', 'k' };

Ii = [ 1 0 0 ; 0 1 0 ; 0 0 1 ; -1 0 0 ; 0 -1 0 ; 0 0 -1 ]';
%Ii = [ 1 0 0 ; 0 1 0 ; 0 0 1 ]';
%Ii = [ 1 0 0 ]';
for k = 1:size(Ii,2)
    for j = 1:3
        I2 = zeros(L, length(range));
        I3 = zeros(3, length(range));

        for i = 1:length(range)
            I2(:,i) = expm(J(:,:,j)*range(i)) * Ii(:,k);

        end

        figure(10+j);
            plot3(I2(1,:), I2(2,:), I2(3,:), sprintf('%c-', char(colors{j})));
            axis([-1 1 -1 1 -1 1]); axis equal;
            hold on;

        figure(4);
            plot3(I2(1,:), I2(2,:), I2(3,:), sprintf('%c-', char(colors{j})));
            axis([-1 1 -1 1 -1 1]); axis equal;
            hold on
    end
end

drawnow;

