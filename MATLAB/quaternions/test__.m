%% PHYSICAL PARAMETERS

d1 = 1; % diameter of particle 1
d2 = 1; % diameter of particle 2
d12 = (d1+d2)/2; % average of diameters

% rho = 1; % density of particles
% m1 = (pi*d1.^3)/6; % mass of particle 1
% I1 = (rho*pi*d1.^5)/40; % moment of inertia of particle 1
% m2 = (pi*d2.^3)/6; % mass of particle 2
% I2 = (rho*pi*d2.^5)/40; % moment of inertia of particle 2

m1 = 1;
m2 = 1;
I1 = (d1.^2)/8;
I2 = (d2.^2)/8;

kd = 1; % constant of dissipaiton
ke = 1; % constant of elasticity

%% INITIAL PARAMETERS

r10 = [0 0 0]; % initial position of particle 1
r20 = [2 0.5 0]; % initial position of particle 2

v10 = [0 0 0]; % initial velocity of particle 1
v20 = [-0.01 0 0]; % initial velocity of particle 2

q10 = [1 0 0 0]; % we chose the initial orientation to be along the x axis
q20 = [1 0 0 0]; % we chose the initial orientation to be along the x axis

w10 = [0 0 0]; % initial rotation vector of particle 1
qp10 = [0 w10(3) -w10(2) -w10(1)]/2; % initial first derivative of first quaternion's derivative
w20 = [0 0 0]; % initial rotation vector of particle 2
qp20 = [0 w20(3) -w20(2) -w20(1)]/2; % initial first derivative of second quaternion's derivative

%% INTEGRATION PARAMETERS

dt = 0.05; % time step
Niter = 4000; % number of iterations so that Niter*dt = time of the simulation

%% INITIALISATION

r1 = [r10;r10 + v10*dt]; % position of particle 1 as a function of time (we assume no forces are exerted at t = 0)
r2 = [r20;r20 + v20*dt]; % position of particle 2 as a function of time (we assume no forces are exerted at t = 0)

v1 = [v10;v10]; % velocity of particle 1 as a funciton of time (we assume no forces are exerted at t = 0)
v2 = [v20;v20]; % velocity of particle 2 as a function of time (we assume no forces are exerted at t = 0)

qp1 = [qp10;qp10 - dt*qnorm(qp10)*q10]; % first derivative of the quaternion of particle 1 as a function of time
q1 = [q10;q10 + dt*qp10]; % quaternion of particle 1 as a function of time, we chose the initial orientation to be along the x axis
qp2 = [qp20;qp20 - dt*qnorm(qp20)*q20]; % first derivative of the quaternion of particle 2 as a function of time
q2 = [q20;q20 + dt*qp20]; % quaternion of particle 2 as a function of time, we chose the initial orientation to be along the x axis

% w1 = [w10;2*quat2vec(qprod(qp1(end,:),qconj(q1(end,:))))];
% w2 = [w20;2*quat2vec(qprod(qp2(end,:),qconj(q2(end,:))))];
w2 = [w20;w20];
w1 = [w10;w10];

%% INTEGRATION

force = [];

for iter = 1:Niter
    
    % useful vectors and scalars
    
    r12 = norm(r1(end,:) - r2(end,:));
    r12vec = r1(end,:) - r2(end,:);
    v12vec = v1(end,:) - v2(end,:);
    s12 = r12vec*(d2/d12);
    s21 = -r12vec*(d1/d12);
    v12 = v1(end,:) - v2(end,:);
    
    % renormalisation
    
    if rem(iter,100) == 0 % the norm of quaternions does not matter, we renormalise them every 100 iterations to avoid erros
        q1(end,:) = q1(end,:)/sqrt(qnorm(q1(end,:)));
        q2(end,:) = q2(end,:)/sqrt(qnorm(q2(end,:)));
    end
    
    % calculation of the forces
    
    if r12 < d12
        fel12 = ke/(d12*r12) * (1 - r12/d12) * r12vec;
        fel21 = - fel12;

        fdis12 = -kd * (v12 + 2*(cross(quat2vec(w2(end,:)),s12) - cross(quat2vec(w1(end,:)),s21)));
        fdis21 = -fdis12;
    else
        fel12 = [0 0 0];
        fel21 = [0 0 0];
        
        fdis12 = [0 0 0];
        fdis21 = [0 0 0];
    end
    force = [force;cat(2,fel12,fdis12)];
    
    % calculation of the moment
    
    M12 = cross(s12,fdis12);
    M21 = cross(s21,fdis21);
    
    % calculation of the accelerations
    
    a1 = (fel12 + fdis12)/m1;
    a2 = (fel21 + fdis21)/m2;
    
    qpp1 = vecscal2quat(qcross(q1(end,:),M12) + q1(end,end)*M12,- dot(M12,quat2vec(q1(end,:))))/(2*I1) - qnorm(qp1(end,:))*q1(end,:);
    qpp2 = vecscal2quat(qcross(q2(end,:),M21) + q2(end,end)*M21,- dot(M21,quat2vec(q2(end,:))))/(2*I2) - qnorm(qp2(end,:))*q2(end,:);
    
    % integration (Verlet's method)
    
    q1 = [q1;2*q1(end,:) - q1(end-1,:) + qpp1*dt.^2];
    q2 = [q2;2*q2(end,:) - q2(end-1,:) + qpp2*dt.^2];
    
    qp1 = [qp1;(q1(end,:) - q1(end-1,:))/dt];
    qp2 = [qp2;(q2(end,:) - q2(end-1,:))/dt];
    
%     w1 = [w1;2*quat2vec(qprod(qp1(end,:),qconj(q1(end,:))))];
%     w2 = [w2;2*quat2vec(qprod(qp2(end,:),qconj(q2(end,:))))];
    w1 = [w1;w1(end,:) + dt*M12/I1];
    w2 = [w2;w2(end,:) + dt*M21/I2];
    
    r1 = [r1;2*r1(end,:) - r1(end-1,:) + a1*dt.^2];
    r2 = [r2;2*r2(end,:) - r2(end-1,:) + a2*dt.^2];
    
    v1 = [v1;(r1(end,:) - r1(end-1,:))/dt];
    v2 = [v2;(r2(end,:) - r2(end-1,:))/dt];
    
end

%% TESTS

r1diff = diff(r1)/dt;
r2diff = diff(r2)/dt;

v1diff = diff(v1)/dt;
v2diff = diff(v2)/dt;

%% DISPLAY

N = 100;

figure;
scatter3(r1(1:N:end,1),r1(1:N:end,2),r1(1:N:end,3));
hold on;
scatter3(r2(1:N:end,1),r2(1:N:end,2),r2(1:N:end,3));

%% SAVE

% Time = 
% x1 = r1(:,1); % position of particle 1 along x axis
% y1 = r1(: