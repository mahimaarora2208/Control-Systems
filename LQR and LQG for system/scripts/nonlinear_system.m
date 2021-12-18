function Z_dot = nonlinear_system(t, y, F, m1, m2, L1, L2, M)
g = 9.81;
x = y(1);
dx = y(2);
t1 = y(3);
dt1 = y(4);
t2 = y(5);
dt2 = y(6);
Z_dot =zeros(6,1);
Z_dot(1) = dx;
Z_dot(2) = (F-((m1*sin(t1)*cos(t1))+(m2*sin(t2)*cos(t2)))*g - (L1*m1*(Z_dot(3)^2)*sin(t1)) - (L2*m2*(Z_dot(5)^2)*sin(t2)))/(m1+m2+M-(m1*(cos(t1)^2))-(m2*(cos(t2)^2)));
Z_dot(3) = dt1;
Z_dot(4) = (cos(t1)*Z_dot(2)-g*sin(t1))/L1;
Z_dot(5) = dt2;
Z_dot(6) = (cos(t2)*Z_dot(2)-g*sin(t2))/L2;
end