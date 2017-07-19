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

w10 = [0 0 0]; % initial rotation vector of particle 1
w20 = [0 0 0]; % initial rotation vector of particle 2

%% INTEGRATION PARAMETERS

dt = 0.05; % time step
Niter = 4000; % number of iterations so that Niter*dt = time of the simulation

%% INITIALISATION

% The integration will resort to Verlet's method of integration, hence
% the 2-step initialisation for positions and velocities and the 1-step
% initialisation for forces.

r1 = [r10;r10 + v10*dt]; % position of particle 1 as a function of time (we assume no force is exerted at t = 0)
r2 = [r20;r20 + v20*dt]; % position of particle 2 as a function of time (we assume no force is exerted at t = 0)

v1 = [v10;v10]; % velocity of particle 1 as a funciton of time (we assume no force is exerted at t = 0)
v2 = [v20;v20]; % velocity of particle 2 as a function of time (we assume no force is exerted at t = 0)

q1 = [0 0 0 1;0 0 0 1]; % quaternion of particle 1 as a function of time, initial q is logically the scalar 1
q2 = [0 0 0 1;0 0 0 1]; % quaternion of particle 2 as a function of time, initial q is logically the scalar 1

qp1 = [vec2quat(w10)/2;vec2quat(w10)/2]; % first derivative of the quaternion of particle 1 as a function of time (we assume no moment is exerted at t = 0)
qp2 = [vec2quat(w20)/2;vec2quat(w20)/2]; % first derivative of the quaternion of particle 2 as a function of time (we assume no moment is exerted at t = 0)

w1 = [w10;w10]; % rotation vector of particle 1 as a function of time (we assume no moment is exerted at t = 0)
w2 = [w20;w20]; % rotation vector of particle 1 as a function of time (we assume no moment is exerted at t = 0)

fel12 = [0 0 0]; % elastic force exerted on particle 1 by partile 2 as a function of time
fel21 = - fel12; % elastic force exerted on particle 2 by partile 1 as a function of time

fdis12 = [0 0 0]; % elastic force exerted on particle 1 by partile 2 as a function of time
fdis21 = - fdis12; % elastic force exerted on particle 2 by partile 1 as a function of time

M12 = [0 0 0]; % moment exerted on particle 1 by particle 2 as a funciton of time
M21 = [0 0 0]; % moment exerted on particle 2 by particle 1 as a function of time

%% INTEGRATION

for iter = 1:Niter

    % useful vectors and scalars
    
    r12 = norm(r1(end,:) - r2(end,:));
    r12vec = r1(end,:) - r2(end,:);
    v12vec = v1(end,:) - v2(end,:);
    s12 = r12vec*(d2/d12);
    s21 = -r12vec*(d1/d12);
    v12 = v1(end,:) - v2(end,:);
    
    % forces at t
    
    if r12 < d12 % particles touch
        
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
    
    % Verlet's method for position
    
    a1 = (fel12(end,:) + fdis12(end,:))/m1; % acceleration felt by particle 1
    a2 = (fel21(end,:) + fdis21(end,:))/m2; % acceleration felt by particle 2

    r1 = [r1;2*r1(end,:) - r1(end-1,:) + a1*dt.^2]; % position of particle 1 at t + dt
    r2 = [r2;2*r2(end,:) - r2(end-1,:) + a2*dt.^2]; % position of particle 2 at t + dt
    
    v1 = [v1;(r1(end,:) - r1(end-1,:))/dt]; % velocity of particle 1 at t + dt
    v2 = [v2;(r2(end,:) - r2(end-1,:))/dt]; % velocity of particle 2 at t + dt
    
    % moments at t
    
    M12 = [M12;cross(s21,fdis12(end,:))]; % moment exerted on particle 1 by particle 2
    M21 = [M21;cross(s12,fdis21(end,:))]; % moment exerted on particle 2 by particle 1
    
    % Kleppmann's method for unitary integration of quaternions
    
    if qp1(end,:) ~= 0
        q1_ = q1(end,:) + tan(dt*sqrt(qnorm(qp1(end,:))))*qp1(end,:)/sqrt(qnorm(qp1(end,:))); % quaternion 1 at t + dt
        q1 = [q1;q1_/sqrt(qnorm(q1_))]; % normalized quaternion 1 at t + dt
    else
        q1 = [q1;q1(end,:)]; % quaternion 1 is unchanged
    end
    
    if qp2(end,:) ~= 0
        q2_ = q2(end,:) + tan(dt*sqrt(qnorm(qp2(end,:))))*qp2(end,:)/sqrt(qnorm(qp2(end,:)));  % quaternion 2 at t + dt
        q2 = [q2;q2_/sqrt(qnorm(q2_))]; % normalized quaternion 2 at t + dt
    else
        q2 =[q2;q2(end,:)]; % quaternion 2 is unchanged
    end
    
    w1 = [w1;w(q1(end,:),qp1(end,:))]; % rotation vector of particle 1
    w2 = [w2;w(q2(end,:),qp2(end,:))]; % rotation vector of particle 2
    
    % Euler's method for the first derivative of quaternion
    
    qp1 = [qp1;qp1(end,:) + dt*qpp(q1(end-1,:),qp1(end,:),M12(end,:),I1)]; % first derivative of the quaternion of particle 1
    qp2 = [qp2;qp2(end,:) + dt*qpp(q2(end-1,:),qp2(end,:),M21(end,:),I2)]; % first derivative of the quaternion of particle 2
    
end

%% DISPLAY

% N = 100;
% 
% figure;
% scatter3(r1(1:N:end,1),r1(1:N:end,2),r1(1:N:end,3));
% hold on;
% scatter3(r2(1:N:end,1),r2(1:N:end,2),r2(1:N:end,3));

%% SAVE

time = transpose(linspace(0,(Niter+1)*dt,Niter+2)); % time corresponding to the values

r10 = transpose(r10); % initial position of particle 1
r20 = transpose(r20); % initial position of particle 2

v10 = transpose(v10); % initial velocity of particle 1
v20 = transpose(v20); % initial velocity of particle 2

w10 = transpose(w10); % initial angular velocity of particle 1
w20 = transpose(w20); % initial angular velocity of particle 2

x1 = r1(:,1); % position of particle 1 along the x axis as a function of time
y1 = r1(:,2); % position of particle 1 along the y axis as a function of time
z1 = r1(:,3); % position of particle 1 along the z axis as a function of time

x2 = r2(:,1); % position of particle 2 along the x axis as a function of time
y2 = r2(:,2); % position of particle 2 along the y axis as a function of time
z2 = r2(:,3); % position of particle 2 along the z axis as a function of time

vx1 = v1(:,1); % velocity of particle 1 along the x axis as a function of time
vy1 = v1(:,2); % velocity of particle 1 along the y axis as a function of time
vz1 = v1(:,3); % velocity of particle 1 along the z axis as a function of time

vx2 = v2(:,1); % velocity of particle 2 along the x axis as a function of time
vy2 = v2(:,2); % velocity of particle 2 along the y axis as a function of time
vz2 = v2(:,3); % velocity of particle 2 along the z axis as a function of time

wx1 = w1(:,1); % angular velocity of particle 1 along the x axis as a function of time
wy1 = w1(:,2); % angular velocity of particle 1 along the y axis as a function of time
wz1 = w1(:,3); % angular velocity of particle 1 along the z axis as a function of time

wx2 = w2(:,1); % angular velocity of particle 2 along the x axis as a function of time
wy2 = w2(:,2); % angular velocity of particle 2 along the y axis as a function of time
wz2 = w2(:,3); % angular velocity of particle 2 along the z axis as a function of time

initial_values = table(r10,r20,v10,v20,w10,w20); % initial values
results = table(time,x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2,wx1,wy1,wz1,wx2,wy2,wz2); % result table

filename = 'results.csv';
writetable(initial_values,filename);
writetable(results,filename);