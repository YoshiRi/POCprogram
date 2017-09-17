% Test for fitting

% define size
N1 = 256;
N2 = N1;
c1 = N1/2;
c2 = N2/2;

% define amplitude and displacement
amp = 1.0;
d1 =12.5;
d2 = 20.2;;
p1 = d1+c1;
p2 = d2+c2;


%% make data
% simpler function
myfunc = @(n1,n2) ( amp/N1/N2 * sin(pi*(n1 + p1)) * sin(pi*(n2 + p2))/sin(pi*(n1+ p1)/N1)/sin(pi*(n2 + p2)/N2));

mat = zeros(N2,N1);
for i = 1:N2
    for j = 1:N1        
        mat(i,j) = myfunc(j,i); 
    end
end

%% show and catch the peak
figure(1);
mesh(mat);

[mm1 x1] = max(mat);
[mx1 y1] = max(mm1);
dx = y1;       % column
dy = x1(y1); % row
display(dx);
display(dy);


%% define Function
syms n1 n2 d1_ d2_ amp_ y

% make sym function
func = amp_/N1/N2 * sin(pi*(n1 + d1_)) * sin(pi*(n2 + d2_))/sin(pi*(n1+ d1_)/N1)/sin(pi*(n2 + d2_)/N2);
R = y - func;

%% crop the map
smat = mat(dy-2:dy+2,dx-2:dx+2);
figure(2);
imshow(smat,[-1 1]);
figure(3);
mesh(smat);

%% fitting on the map
%
data_y  = [1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5].' + dy-3;
data_x  = [1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5].' + dx-3;
data = data_x;
for i = 1:length(data_y)
    data(i) = smat(i);
end

% partial derivative
fd1 = diff(R,d1_);
fd2 = diff(R,d2_);
famp = diff(R,amp_);


% make jacobian matrix
Jd1 = subs(fd1,{n1,n2},{data_x,data_y});
Jd2 = subs(fd2,{n1,n2},{data_x,data_y});
Jamp = subs(famp,{n1,n2},{data_x,data_y});

Jacob = [Jd1,Jd2,Jamp]; % jacobian
Jacob2 =[Jd1,Jd2];
r = subs(R,{n1,n2,y},{data_x,data_y,data});
r2 = subs(R,{n1,n2,y,amp_},{data_x,data_y,data,ones(size(data))});


d1_hat = N1- dx +1;
d2_hat = N2 - dy +1;
amp_hat = mx1;
% kousin
v = [d1_hat, d2_hat, amp_hat].';
v2 = [d1_hat, d2_hat].';

% lambda
lambda = 1;
I = eye(size(Jacob,2));
I2 = eye(size(Jacob2,2));

for i = 1:10
J = double(subs(Jacob,{d1_,d2_,amp_},{v(1),v(2),v(3)}));
J2 = double(subs(Jacob2,{d1_,d2_,amp_},{v2(1),v2(2),amp}));

% differ
df = double(subs(r,{d1_,d2_,amp_},{v(1),v(2),v(3)}));
df2 = double(subs(r2,{d1_,d2_},{v2(1),v2(2)}));

% 更新量Δを計算して更新
delta = - ( J.' * J )\J.' *  df;
delta2 = - ( J2.' * J2 + lambda*I2 )\J2.' *  df2;


dmove = delta./abs(delta)/40;
dmove2 = delta2./abs(delta2)/40;

error = df.' * df
error2 = df2.' * df2

v = v+dmove;
v2 = v2+dmove2;

end



%% Mapping
sX = dx-1:0.1:dx+1;
sY = dy-1:0.1:dy+1;
len = length(sX);

Map = zeros(len);

for i = 1:len
    for j = 1:len
       Map(i,j) = myfunc(sX(j),sY(i)); 
    end
end

figure(3);
mesh(Map);


%% gradMap

gMap1 = Map;
gMap2 = Map;

for i = 1:len
    for j = 1:len
       gMap1(i,j) = subs(fd1,[amp_,n1,n2,d1_,d2_],[1,sX(j),sY(i),dx,dy]); 
       gMap2(i,j) = subs(fd2,[amp_,n1,n2,d1_,d2_],[1,sX(j),sY(i),dx,dy]); 
    end
end

