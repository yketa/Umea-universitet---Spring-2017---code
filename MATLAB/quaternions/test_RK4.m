%% PHYSICAL PARAMETERS

d1 = 1; % diameter of particle 1
d2 = 1; % diameter of particle 2
d12 = (d1+d2)/2; % average of diameters

m1 = 1; % mass of particle 1
m2 = 1; % mass of particle 2
I1 = (d1.^2)/8; % moment of inertia of particle 1
I2 = (d2.^2)/8; % moment of inertia of particle 2

kd = 1; % constant of dissipaiton
ke = 1; % constant of elasticity

%% INITIAL PARAMETERS

% We will consider that are particles are not rotating at t = 0.

r10 = [0 0 0]; % initial position of particle 1
r20 = [2 0.5 0]; % initial position of particle 2

v10 = [0 0 0]; % initial velocity of particle 1
v20 = [-0.01 0 0]; % initial velocity of particle 2

%% INTEGRATION PARAMETERS

dt = 0.05; % time step
Niter = 4000; % number of iterations so that Niter*dt = time of the simulation

%% INITIALISATION

% The integration will resort to cubic Lagrange polynomial with 3 steps of
% history, hence the 3-step initialisation for positions and velocities and
% the 2-step initialisation for forces.

r1 = [r10;r10 + v10*dt;r10 + 2*v10*dt]; % position of particle 1 as a function of time (we assume no force is exerted at t = 0)
r2 = [r20;r20 + v20*dt;r20 + 2*v20*dt]; % position of particle 2 as a function of time (we assume no force is exerted at t = 0)

v1 = [v10;v10;v10]; % velocity of particle 1 as a funciton of time (we assume no force is exerted at t = 0)
v2 = [v20;v20;v20]; % velocity of particle 2 as a function of time (we assume no force is exerted at t = 0)

q1 = [0 0 0 1;0 0 0 1;0 0 0 1]; % quaternion of particle 1 as a function of time, initial q is logically the scalar 1
q2 = [0 0 0 1;0 0 0 1;0 0 0 1]; % quaternion of particle 2 as a function of time, initial q is logically the scalar 1

w1 = [0 0 0;0 0 0;0 0 0]; % rotation vector of particle 1 as a funciton of time (we assume no moment is exerted at t = 0)
w2 = [0 0 0;0 0 0;0 0 0]; % rotation vector of particle 2 as a funciton of time (we assume no moment is exerted at t = 0)

qp1 = [0 0 0 0;0 0 0 0;0 0 0 0]; % first derivative of the quaternion of particle 1 as a function of time
qp2 = [0 0 0 0;0 0 0 0;0 0 0 0]; % first derivative of the quaternion of particle 1 as a function of time

fel12 = [0 0 0;0 0 0]; % elastic force exerted on particle 1 by partile 2 as a function of time
fel21 = - fel12; % elastic force exerted on particle 2 by partile 1 as a function of time

fdis12 = [0 0 0;0 0 0]; % elastic force exerted on particle 1 by partile 2 as a function of time
fdis21 = - fdis12; % elastic force exerted on particle 2 by partile 1 as a function of time

M12 = [0 0 0;0 0 0]; % moment exerted on particle 1 by particle 2 as a funciton of time
M21 = [0 0 0;0 0 0]; % moment exerted on particle 2 by particle 1 as a function of time

%% INTEGRATION

prout =[];

for iter = 1:Niter
        
    % renormalisation
    
%     q1(end,:) = q1(end,:)/sqrt(qnorm(q1(end,:)));
%     q2(end,:) = q2(end,:)/sqrt(qnorm(q2(end,:)));
    
%     if rem(iter,100) == 0 % the norm of quaternions does not matter, we renormalise them every 100 iterations to avoid erros
%         q1(end,:) = q1(end,:)/s
%     end

    % useful vectors and scalars
    
    r12 = norm(r1(end,:) - r2(end,:));
    r12vec = r1(end,:) - r2(end,:);
    v12vec = v1(end,:) - v2(end,:);
    s12 = r12vec*(d2/d12);
    s21 = -r12vec*(d1/d12);
    v12 = v1(end,:) - v2(end,:);
    
    % forces at t
    
    if r12 < d12 % particles are in contact
        
        fel12 = [fel12;ke/(d12*r12) * (1 - r12/d12) * r12vec]; % elastic force exerted on particle 1 by particle 2
        fel21 = [fel21;- fel12(end,:)]; % elastic force exerted on particle 2 by particle 1

        fdis12 = [fdis12;-kd * (v12 + cross(quat2vec(w1(end,:)),s21) - cross(quat2vec(w2(end,:)),s12))]; % dissipation force exerted on particle 1 by particle 2
        fdis21 = [fdis21;- fdis12(end,:)]; % dissipation force exerted on particle 2 by particle 1
        
    else % particles do not touch
        
        fel12 = [fel12;0 0 0]; % elastic force exerted on particle 1 by particle 2
        fel21 = [fel21;0 0 0]; % elastic force exerted on particle 2 by particle 1
        
        fdis12 = [fdis12;0 0 0]; % dissipation force exerted on particle 1 by particle 2
        fdis21 = [fdis21;0 0 0]; % dissipation force exerted on particle 2 by particle 1
        
    end
     
    % moments at t
    
    M12 = [M12;cross(s21,fdis12(end,:))]; % moment exerted on particle 1 by particle 2
    M21 = [M21;cross(s12,fdis21(end,:))]; % moment exerted on particle 2 by particle 1
    
    % extrapolated forces at t' > t (cubic Lagrange polynomial)
    
    fel12_ = (3*fel12(end-2,:))/8 - (5*fel12(end-1,:))/4 + (15*fel12(end,:))/8; % extrapolation at t + dt/2
    fel12__ = fel12(end-2,:) - 3*fel12(end-1,:) + 3*fel12(end,:); % extrapolation at t + dt
    
    fel21_ = - fel12_; % etrapolation at t + dt/2
    fel21__ = - fel12__; % extrapolation at t + dt
    
    fdis12_ = (3*fdis12(end-2,:))/8 - (5*fdis12(end-1,:))/4 + (15*fdis12(end,:))/8; % extrapolation at t + dt/2
    fdis12__ = fdis12(end-2,:) - 3*fdis12(end-1,:) + 3*fdis12(end,:); % extrapolation at t + dt
    
    fdis21_ = - fdis12_; % etrapolation at t + dt/2
    fdis21__ = - fdis12__; % extrapolation at t + dt
    
    % extrapolated moments at t' > t (cubic Lagrange polynomial)
    
    M12_ = (3*M12(end-2,:))/8 - (5*M12(end-1,:))/4 + (15*M12(end,:))/8; % extrapolation at t + dt/2
    M12__ = M12(end-2,:) - 3*M12(end-1,:) + 3*M12(end,:); % extrapolation at t + dt
    prout = [prout;cat(2,M12_,M12__)];
    
    M21_ = (3*M21(end-2,:))/8 - (5*M21(end-1,:))/4 + (15*M21(end,:))/8; % extrapolation at t + dt/2
    M21__ = M21(end-2,:) - 3*M21(end-1,:) + 3*M21(end,:); % extrapolation at t + dt

    % 1st Runge-Kutta increments for quaternions
    
    k1q1 = qpp(q1(end,:),qp1(end,:),M12(end,:),I1);
    k1q2 = qpp(q2(end,:),qp2(end,:),M21(end,:),I2);
    
    k2q1 = qpp(q1(end,:) + (dt/2)*qp1(end,:) + (dt.^2/8)*k1q1,qp1(end,:) + (dt/2)*k1q1,M12_,I1);
    k2q2 = qpp(q2(end,:) + (dt/2)*qp2(end,:) + (dt.^2/8)*k1q2,qp2(end,:) + (dt/2)*k1q2,M21_,I2);
    
    k3q1 = qpp(q1(end,:) + (dt/2)*qp1(end,:) + (dt.^2/8)*k1q1,qp1(end,:) + (dt/2)*k2q1,M12_,I1);
    k3q2 = qpp(q2(end,:) + (dt/2)*qp2(end,:) + (dt.^2/8)*k1q2,qp2(end,:) + (dt/2)*k2q2,M21_,I2);
    
    k4q1 = qpp(q1(end,:) + dt*qp1(end,:) + (dt.^2/2)*k3q1,qp1(end,:) + dt*k3q1,M12__,I1);
    k4q2 = qpp(q2(end,:) + dt*qp2(end,:) + (dt.^2/2)*k3q2,qp2(end,:) + dt*k3q2,M21__,I2);
    
    % 1st Runge-Kutta increments for positions
    
    k1r1 = (fel12(end,:) + fdis12(end,:))/m1;
    k1r2 = (fel21(end,:) + fdis21(end,:))/m2;
    
    k2r1 = (fel12_ + fdis12_)/m1;
    k2r2 = (fel21_ + fdis21_)/m2;
    
    k3r1 = (fel12__ + fel12__)/m1;
    k3r2 = (fel21__ + fel21__)/m2;
    
    % Kleppmann's method for unitary integration of quaternions
    
    if qp1(end,:) ~= 0
        q1_ = q1(end,:) + tan(dt*sqrt(qnorm(qp1(end,:))))*qp1(end,:)/sqrt(qnorm(qp1(end,:)));
        q1 = [q1;q1_/sqrt(qnorm(q1_))];
    else
        q1 = [q1;q1(end,:)];
    end
    
    if qp2(end,:) ~= 0
        q2_ = q2(end,:) + tan(dt*sqrt(qnorm(qp2(end,:))))*qp2(end,:)/sqrt(qnorm(qp2(end,:)));
        q2 = [q2;q2_/sqrt(qnorm(q2_))];
    else
        q2 =[q2;q2(end,:)];
    end
    
    w1 = [w1;w(q1(end,:),qp1(end,:))];
    w2 = [w2;w(q2(end,:),qp2(end,:))]; 
    
    % Runge-Kutta integration for the first derivatives of quaternions
    
%     q1 = [q1;q1(end,:) + dt*qp1(end,:) + (dt.^2/6)*(k1q1 + k2q1 + k3q1)];
%     q2 = [q2;q2(end,:) + dt*qp2(end,:) + (dt.^2/6)*(k1q2 + k2q2 + k3q2)];
    
    qp1 = [qp1;qp1(end,:) + (dt/6)*(k1q1 + 2*k2q1 + 2*k3q1 + k4q1)];
    qp2 = [qp2;qp2(end,:) + (dt/6)*(k1q2 + 2*k2q2 + 2*k3q2 + k4q2)];
    
%     w1 = [w1;2*quat2vec(qprod(qp1(end,:),qconj(q1(end,:))))];
%     w2 = [w2;2*quat2vec(qprod(qp2(end,:),qconj(q2(end,:))))];

    % Runge-Kutta integration for position
    
    r1 = [r1;r1(end,:) + dt*v1(end,:) + (dt.^2/6)*(k1r1 + 2*k2r1)];
    r2 = [r2;r2(end,:) + dt*v2(end,:) + (dt.^2/6)*(k1r2 + 2*k2r2)];
    
    v1 = [v1;v1(end,:) + (dt/6)*(k1r1 + 4*k2r1 + k3r1)];
    v2 = [v2;v2(end,:) + (dt/6)*(k1r2 + 4*k2r2 + k3r2)];
    
end

%% TESTS

r1diff = diff(r1)/dt;
r2diff = diff(r2)/dt;

v1diff = diff(v1)/dt;
v2diff = diff(v2)/dt;

qn1 = qnorm(q1);
qn2 = qnorm(q2);

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