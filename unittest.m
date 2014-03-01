
%% I1 = expm(J*t) * I0
%% E = ||I1 - expm(J*t)*I0||_2^2 + lambda*||t||_2^2 + gamma*||J||_2^2

clear

L = 3;
M = 9;
Msz = sqrt(M);

display_every = 100;
save_every = 100;

shift_rad = 0.5;
flag_display = 1;
flag_learn = 1;

lambda = 0.0001;  % l1 prior on c
gamma = 0.01;    % l2 prior on J
dangle = 0.01;

eta_J = 0.05;


paramstr = sprintf('L=%03d_%s',L,datestr(now,30));

reinit

batch_size = 1;
num_trials = 1000000;
mapnet


