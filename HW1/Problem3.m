clc; clear; close all;

%time as a variable
syms t;

%Gravitational constant
syms g;
g=9.81;

%masses (m1 is unused, but placed for future potential use in different
%configurations)
syms m1 m2 m3;

%Dimensional constants
syms a1 a2 a3;
syms l1 l2 l3
l1 = a1/2
l2 = a2/2
l3 = a3/3

%Rotational inertias, found by parallel axis theorum
syms I1 I2 I3
I1 = m1*l1^2
I2 = m2*l2^2
I3 = m3*l3^3

%Joint variables. Rotational displacement in radians
syms q1(t) q2(t) q3(t);

%Position and velocity vectors of center of masses of links
position_1 = [l1*sin(q1);l1*cos(q1)]
position_2 = [a1*sin(q1)+l2*sin(q1+q2);a1*cos(q1)+l2*cos(q1+q2)]
position_3 = [a1*sin(q1)+a2*sin(q1+q2)+l3*sin(q1+q2-q3);
              a1*cos(q1)+a2*cos(q1+q2)+l3*cos(q1+q2-q3)]

%Used for indexing symbolic matrices later in code
position1symbolic = position_1(t)
position2symbolic = position_2(t)
position3symbolic = position_3(t)

velocity_1 = diff(position_1,t)
velocity_2 = diff(position_2,t)
velocity_3 = diff(position_3,t)

v1squared = transpose(velocity_1)*velocity_1
v2squared = transpose(velocity_2)*velocity_2
v3squared = transpose(velocity_3)*velocity_3

%Used for indexing symbolic matrices later in code
v1squaredsymbolic = v1squared(t)
v2squaredsymbolic = v2squared(t)
v3squaredsymbolic = v3squared(t)

%Kinetic energy contributed from each joint is equal to the linear velocity
%of the center of mass, plus the rotational inertia of the link with its 
%total angular velocity from the joint and all prior joints

Ke1 = 0.5*m1*v1squaredsymbolic+0.5*I1*(diff(q1,t))^2
Ke2 = 0.5*m2*v2squaredsymbolic+0.5*I2*(diff(q2,t))^2
Ke3 = 0.5*m3*v3squaredsymbolic+0.5*I3*(diff(q3,t))^2

%Potential energy from link height above Y=0
U1 = m1*g*position1symbolic(2)
U2 = m2*g*position1symbolic(2)
U3 = m3*g*position1symbolic(2)

%Lagrangian, total of K-U
L=Ke1+Ke2+Ke3-U1-U2-U3;

%To get derivatives with respect to joint velocity, diff(qi(t),t) must be 
%replaced with qdot variables
syms q1dot(t) q2dot(t) q3dot(t);
Lsymbolic1 = subs(L,diff(q1(t), t),q1dot);
Lsymbolic2 = subs(Lsymbolic1,diff(q2(t), t),q2dot);
Lsymbolic3 = subs(Lsymbolic2,diff(q3(t), t),q3dot);
L = Lsymbolic3

%Solve for each joint torque using d/dt dL/dqdoti - dL/dqi
syms T1 T2 T3;
%Find derivative with repsect to joint velocity
dqdot1 = functionalDerivative(L, q1dot(t));
dqdot2 = functionalDerivative(L, q2dot(t));
dqdot3 = functionalDerivative(L, q3dot(t));
T1 = diff(dqdot1, t) - functionalDerivative(L, q1)
T2 = diff(dqdot2, t) - functionalDerivative(L, q2)
T3 = diff(dqdot3, t) - functionalDerivative(L, q3)









