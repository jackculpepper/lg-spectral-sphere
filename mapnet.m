

for t=1:num_trials

    tic


    %% generate points on a sphere
    I0 = zeros(L, batch_size);
    I1 = zeros(L, batch_size);

    psi0 = zeros(batch_size);
    theta0 = zeros(batch_size);
    psi1 = zeros(batch_size);
    theta1 = zeros(batch_size);

    for i = 1:batch_size
        theta0(i) = pi*rand;
        psi0(i) = 2*pi*rand;

        I0(:,i) = [ sin(theta0(i))*cos(psi0(i)) ; ...
                    sin(theta0(i))*sin(psi0(i)) ; ...
                    cos(theta0(i)) ];

        theta1(i) = theta0(i) + dangle*(randn);
        psi1(i) = psi0(i) + dangle*(randn);

        I1(:,i) = [ sin(theta1(i))*cos(psi1(i)) ; ...
                    sin(theta1(i))*sin(psi1(i)) ; ...
                    cos(theta1(i)) ];
    end


    time_choose = toc;

    tic

    %% init solution
    c0 = 0.1*rand(M, batch_size);
    c0 = zeros(M, batch_size);
    c1 = zeros(M, batch_size);

    %% run inference
    for i = 1:batch_size
        if 0
            figure(5);
            subplot(2,2,1); imagesc(reshape(I0, Lsz, Lsz)); colormap(gray);
                axis image; colorbar;
            subplot(2,2,2); imagesc(reshape(I1, Lsz, Lsz)); colormap(gray);
                axis image; colorbar;
                title(sprintf('dx %.4f dy %.4f', size(F0,1)*dx/shift_rad, ...
                                                 size(F0,2)*dy/shift_rad));
        end

        c1(:,i) = minimize(c0(:,i), 'objfun_c_spectral_fast', 20, J, ...
                           I0(:,i), I1(:,i), lambda);
    end
    c1 = real(c1);

    if (batch_size == 1)
        figure(6);
        bar(c1); title('c');
    end

    time_inf = toc;

    %% generate images from solution
    EI1 = zeros(L, batch_size);
    A = zeros(L,L);
    for i = 1:M
        A = A + J(:,:,i) * c1(i);
    end
    ExpA = expm(A);
    for i = 1:batch_size
        EI1(:,i) = ExpA * I0(:,i);
    end

    %% compare actual and estimated ; compute snr
    E = I1 - EI1;
    snr = 10 * log10 ( sum(I1(:).^2) / sum(E(:).^2) );

    if 0
        figure(5);
        subplot(2,2,3); imagesc(reshape(ExpA*I0, Lsz, Lsz)); colormap(gray);
            axis image; colorbar;
        subplot(2,2,4); imagesc(reshape(E, Lsz, Lsz)); colormap(gray);
            axis image; title(sprintf('snr %.4f', snr)); colorbar;
        drawnow;
    end
    
    tic


    %% accumulate fJ and dJ
    fJ0 = 0;
    J1 = J;
    for i = 1:batch_size
        [fJ, dJ] = objfun_J_mult(J(:),c1(:,i),I0(:,i),I1(:,i),gamma);

        fJ0 = fJ0 + real(fJ);
        J1 = J1 - eta_J * reshape(real(dJ), L, L, M);
    end

    %% accumulate fJ after the update
    fJ1 = 0;
    for i = 1:batch_size
        [fJ, dJ] = objfun_J_mult(J1(:),c1(:,i),I0(:,i),I1(:,i),gamma);
        fJ1 = fJ1 + real(fJ);
    end

    if flag_learn

        %% take the update if fJ went down
        if fJ1 < fJ0
            J = J1;
        else
            fprintf('warning: objfun_J increased: skipping update\n');
            eta_J = eta_J * 0.99;
        end
    end

    time_updt = toc;

    % display
    
    tic
    if (display_every == 1 || mod(update,display_every)==0)
        fprintf('rendering display..');

        figure(1); colormap(gray);
        for i = 1:M
            subplot(Msz,Msz,i); imagesc(J(:,:,i)); axis image; colorbar;
        end

        figure(3);
        Ii = [1 0 0]';
        range = -1000:10:1000;
        for j = 1:M
            I2 = zeros(L, length(range));
            for i = 1:length(range)
                I2(:,i) = expm(J(:,:,j)*range(i)) * Ii;
            end

            subplot(Msz,Msz,j);

            plot3(I2(1,:), I2(2,:), I2(3,:), '.');
            axis([-1 1 -1 1 -1 1]);
        end

        drawnow;

        if (save_every == 1 || mod(update,save_every)==0)

            [sucess,msg,msgid]=mkdir(sprintf('state/%s', paramstr));
 
            if 0
                array_frame = uint8(255*((array+1)/2)+1);
                imwrite(array_frame, ...
                    sprintf('state/%s/l1bf_up=%06d.gif', ...
                    paramstr,update), 'gif');
            end

            eval(sprintf('save state/%s/J.mat J', paramstr));

            saveparamscmd = sprintf('save state/%s/params.mat', paramstr);
            saveparamscmd = sprintf('%s lambda', saveparamscmd);
            saveparamscmd = sprintf('%s gamma', saveparamscmd);
            saveparamscmd = sprintf('%s eta_J', saveparamscmd);
            saveparamscmd = sprintf('%s L', saveparamscmd);
            saveparamscmd = sprintf('%s M', saveparamscmd);
            eval(saveparamscmd);

        end
        fprintf('done\r');
    end
    time_disp = toc;


    fprintf('update %d', update);
    fprintf(' %s', paramstr);
    %fprintf(' dx %.2f dy %.2f', dx, dy);
    fprintf(' fJ0 %.2f fJ1 %.2f', fJ0, fJ1);
    fprintf(' snr %.4f\n', snr);

    %% renormalize J
    %for i = 1:M
    %    J(:,:,i) = J(:,:,i)*diag(1./sqrt(sum(sum(J(:,:,i).^2))));
    %end



    update = update + 1;
end


