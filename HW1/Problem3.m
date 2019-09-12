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

%Joint variables. Rotational displacement in radians
syms q1(t) q2(t) q3(t);

%Position and velocity vectors
position_1 = [a1*sin(q1);a1*cos(q1)]
position_2 = [a1*sin(q1)+a2*sin(q1+q2);a1*cos(q1)+a2*cos(q1+q2)]
position_3 = [a1*sin(q1)+a2*sin(q1+q2)+a3*sin(q1+q2-q3);
              a1*cos(q1)+a2*cos(q1+q2)+a3*cos(q1+q2-q3)]

%Velocity squared, for kinetic energy
vsquared = transpose(velocity_t)*velocity_t;

%Rotational inertia will be found via the parallel axis theorum, based on
%the distance between mt and the z axis
syms I
I = mt*((b+q2)^2);

%Kinetic energy is the sum of the linear energy in each direction, plus the
%rotational energy from the rotation of q1
Ke = 0.5*mt*vsquared + 0.5*I*(diff(q1,t))^2;

%Potential energy measured as height from Z=0
possymbolic=position_t(t);
U = mt*g*possymbolic(3)

%Compute lagrangian
L=Ke-U


%To get derivatives with respect to joint velocity, diff(qi(t),t) must be 
%replaced with qdot variables
syms q1dot(t) q2dot(t) q3dot(t);
Lsymbolic1 = subs(L,diff(q1(t), t),q1dot);
Lsymbolic2 = subs(Lsymbolic1,diff(q2(t), t),q2dot);
Lsymbolic3 = subs(Lsymbolic2,diff(q3(t), t),q3dot);
L = Lsymbolic3;

%Solve for each joint torque using d/dt dL/dqdoti - dL/dqi
syms T1 T2 T3;
%Find derivative with repsect to joint velocity
dqdot1 = functionalDerivative(L, q1dot(t));
dqdot2 = functionalDerivative(L, q2dot(t));
dqdot3 = functionalDerivative(L, q3dot(t));
T1 = diff(dqdot1, t) - functionalDerivative(L, q1)
T2 = diff(dqdot2, t) - functionalDerivative(L, q2)
T3 = diff(dqdot3, t) - functionalDerivative(L, q3)


