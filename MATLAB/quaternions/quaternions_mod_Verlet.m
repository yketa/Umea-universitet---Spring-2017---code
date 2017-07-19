%% PHYSICAL PARAMETERS

d1 = 1; % diameter of particle 1
d2 = 1; % diameter of particle 2
d12 = (d1+d2)/2; % average of diameters

m1 = 1; % mass of particle 1
m2 = 1; % mass of particle 2

I1 = m1*(d1.^2)/8; % moment of inertia of particle 1
I2 = m2*(d2.^2)/8; % moment of inertia of particle 2

kd = 1; % constant of dissipaiton
ke = 1; % constant of elasticity

%% INITIAL PARAMETERS

% We will consider that are particles are not rotating at t = 0.

r10 = [0 0 0]; % initial position of particle 1
r20 = [2 0.5 0.5]; % initial position of particle 2

v10 = [0 0 0]; % initial velocity of particle 1
v20 = [-0.01 0 0]; % initial velocity of particle 2

w10 = [1 0 0]; % initial rotation vector of particle 1
w20 = [0 0 0]; % initial rotation vector of particle 2

%% INTEGRATION PARAMETERS

dt = 0.05; % time step
Niter = 4000; % number of iterations so that Niter*dt = time of the simulation

%% INITIALISATION

r1 = [r10]; % position of particle 1 as a function of time
r2 = [r20]; % position of particle 2 as a function of time

v1 = [v10]; % velocity of particle 1 as a funciton of time
v2 = [v20]; % velocity of particle 2 as a function of time

q1 = [0 0 0 1]; % quaternion of particle 1 as a function of time, initial q is logically the scalar 1
q2 = [0 0 0 1]; % quaternion of particle 2 as a function of time, initial q is logically the scalar 1

qp1 = [vec2quat(w10)/2]; % first derivative of the quaternion of particle 1 as a function of time
qp2 = [vec2quat(w20)/2]; % first derivative of the quaternion of particle 2 as a function of time

w1 = [w10]; % rotation vector of particle 1 as a function of time
w2 = [w20]; % rotation vector of particle 1 as a function of time

%% INTEGRATION

for iter = 1:Niter

    % forces at t
    
    f1 = fel(ke,r1(end,:),d1,r2(end,:),d2) + fdis(kd,r1(end,:),v1(end,:),w1(end,:),d1,r2(end,:),v2(end,:),w2(end,:),d2); % net force exerted on particle 1 by particle 2
    f2 = fel(ke,r2(end,:),d2,r1(end,:),d1) + fdis(kd,r2(end,:),v2(end,:),w2(end,:),d2,r1(end,:),v1(end,:),w1(end,:),d1); % net force exerted on particle 2 by particle 1
    
    % moments at t
    
    M1 = Mdis(kd,r1(end,:),v1(end,:),w1(end,:),d1,r2(end,:),v2(end,:),w2(end,:),d2); % moment exerted on particle 1 by particle 2
    M2 = Mdis(kd,r2(end,:),v2(end,:),w2(end,:),d2,r1(end,:),v1(end,:),w1(end,:),d1); % moment exerted on particle 1 by particle 2
    
    % first increments
    
    vt11 = v1(end,:) + f1*(dt/m1); % first increment in velocity for particle 1
    vt21 = v2(end,:) + f2*(dt/m2); % first increment in velocity for particle 2
    
    rt11 = r1(end,:) + v1(end,:)*dt; % first increment in position for particle 1
    rt21 = r2(end,:) + v2(end,:)*dt; % first increment in position for particle 2
    
    qpt11 = qp1(end,:) + qpp(q1(end,:),qp1(end,:),M1,I1)*dt; % first increment in first derivative of the quaternion for particle 1
    qpt21 = qp2(end,:) + qpp(q2(end,:),qp2(end,:),M2,I2)*dt; % first increment in first derivative of the quaternion for particle 2
    
    qt11 = qunit(q1(end,:) + tan(dt*sqrt(qnorm(qp1(end,:))))*qunit(qp1(end,:))); % first increment in quaternion for particle 1
    qt21 = qunit(q2(end,:) + tan(dt*sqrt(qnorm(qp2(end,:))))*qunit(qp2(end,:))); % first increment in quaternion for particle 2
    
    wt11 = w(qt11,qpt11); % first increment in angular velocity for particle 1
    wt21 = w(qt21,qpt21); % first increment in angular velocity for particle 2
    
    % forces at first increment
    
    f1_ = fel(ke,rt11,d1,rt21,d2) + fdis(kd,rt11,vt11,wt11,d1,rt21,vt21,wt21,d2); % net force exerted on particle 1 by particle 2 at first increment
    f2_ = fel(ke,rt21,d2,rt11,d1) + fdis(kd,rt21,vt21,wt21,d2,rt11,vt11,wt11,d1); % net force exerted on particle 2 by particle 1 at first increment
    
    % moments at first increment
    
    M1_ = Mdis(kd,rt11,vt11,wt11,d1,rt21,vt21,wt21,d2); % moment exerted on particle 1 by particle 2 at first increment
    M2_ = Mdis(kd,rt21,vt21,wt21,d2,rt11,vt11,wt11,d1); % moment exerted on particle 1 by particle 2 at first increment
    
    % integration
    
    v1 = [v1;v1(end,:) + (f1 + f1_)*(dt/(2*m1))]; % velocity of particle 1 at t + dt
    v2 = [v2;v2(end,:) + (f2 + f2_)*(dt/(2*m2))]; % velocity of particle 2 at t + dt
    
    r1 = [r1;r1(end,:) + (vt11 + v1(end,:))*(dt/2)]; % position of particle 1 at t + dt
    r2 = [r2;r2(end,:) + (vt21 + v2(end,:))*(dt/2)]; % position of particle 1 at t + dt
    
    qp1 = [qp1;qp1(end,:) + (qpp(q1(end,:),qp1(end,:),M1,I1) + qpp(qt11,qpt11,M1_,I1))*(dt/2)]; % first derivative of the quaternion of particle 1 at t + dt
    qp2 = [qp2;qp2(end,:) + (qpp(q2(end,:),qp2(end,:),M2,I2) + qpp(qt21,qpt21,M2_,I2))*(dt/2)]; % first derivative of the quaternion of particle 2 at t + dt
    
    q1inc = (qp1(end,:) + qpt11)*(dt/2); % increment in quaternion for particle 1 at t
    q2inc = (qp2(end,:) + qpt21)*(dt/2); % increment in quaternion for particle 2 at t
    
    q1 = [q1;qunit(q1(end,:) + tan(sqrt(qnorm(q1inc)))*qunit(q1inc))]; % quaternion of particle 1 at t + dt
    q2 = [q2;qunit(q2(end,:) + tan(sqrt(qnorm(q2inc)))*qunit(q2inc))]; % quaternion of particle 1 at t + dt
    
    w1 = [w1;w(q1(end,:),qp1(end,:))]; % rotation vector of particle 1 at t + dt
    w2 = [w2;w(q2(end,:),qp2(end,:))]; % rotation vector of particle 1 at t + dt
    
end

%% DISPLAY

N = 100;

figure;
scatter3(r1(1:N:end,1),r1(1:N:end,2),r1(1:N:end,3));
hold on;
scatter3(r2(1:N:end,1),r2(1:N:end,2),r2(1:N:end,3));

%% SAVE

time = transpose(linspace(0,Niter*dt,Niter+1)); % time corresponding to the values

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

Etrans1 = Etrans(v1,m1); % translationnal energy of particle 1
Etrans2 = Etrans(v2,m2); % translationnal energy of particle 2
Etranstot = Etrans1 + Etrans2; % total translationnal energy

Erot1 = Erot(w1,I1); % rotationnal energy of particle 1
Erot2 = Erot(w2,I2); % rotationnal energy of particle 2
Erottot = Erot1 + Erot2; % total rotationnal energy

Etot = Etranstot + Erottot; % total energy

results = table(time,x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2,wx1,wy1,wz1,wx2,wy2,wz2,Etrans1,Etrans2,Etranstot,Erot1,Erot2,Erottot,Etot); % result table

filename = 'results.csv';
writetable(results,filename);