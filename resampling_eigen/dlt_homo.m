h11=sym('h11', 'real');
h12=sym('h12', 'real');
h13=sym('h13', 'real');
h21=sym('h21', 'real');
h22=sym('h22', 'real');
h23=sym('h23', 'real');
h31=sym('h31', 'real');
h32=sym('h32', 'real');
h33=sym('h33', 'real');
x1=sym('x1', 'real');
x2=sym('x2', 'real');
w = sym('w', 'real');

y1=sym('y1', 'real');
y2=sym('y2', 'real');
xk = [(w*x1) ; (w*x2) ; w]
yk = [y1 ; y2; 1]

H = [h11 h12 h13 ; h21 h22 h23 ; h31 h32 h33]

w_tilde = H(3, 1:3) * yk
wx1_tilde = H(1, 1:3)*yk
wx2_tilde = H(2, 1:3)*yk

zx1 = (w_tilde * x1) -  wx1_tilde
zx2 = (w_tilde * x2) -  wx2_tilde


%Mise en forme de A pour consituer la matrice ...

Ak =  [-yk(1), -yk(2), -1, 0, 0, 0, yk(1)*x1, yk(2)*x1 , x1;
        0  ,  0,0, -yk(1),-yk(2),-1, yk(1)*x2, yk(2)*x2, x2 ] 

%Must concat Ak for k -> 1 to 4 and then use SVD

%for each point pair
% the solution of the system is th


%now need to create a set of 4 points
% apply homography to these points
% compute homography with aforementionned method

