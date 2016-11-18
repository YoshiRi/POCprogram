function [amp pxx pyy] = POCfitting(peak_x ,peak_y, Pp)

% initial value
d1 = peak_y;
d2 = peak_x;
alpha = Pp(d1,d2);

% parameter
N1 = size(Pp,2);
N2 = size(Pp,1);

% make refer data
y = [Pp(d1,d2) Pp(d1-1,d2) Pp(d1+1,d2) Pp(d1,d2+1) Pp(d1-1,d2+1) Pp(d1+1,d2+1) Pp(d1,d2-1) Pp(d1-1,d2-1) Pp(d1+1,d2-1)].';
x1 = [d1 d1+1 d1-1 d1 d1+1 d1-1 d1 d1+1 d1-1 ].'; 
x2 = [d2 d2 d2 d2+1 d2+1 d2+1 d2-1 d2-1 d2-1].';

% make initial lambda
lambda = init_lambda(y,x1,x2,d1,d2,alpha,N1,N2);
init_cost = cost_function_POC(y,x1,x2,d1,d2,alpha,N1,N2);
p_cost = init_cost;

% begin roop
maxroop = 25;
finish_magnitude = 1/10;

while maxroop > 0
    delta = leastsquare_Marquardt(y,x1,x2,d1,d2,alpha,N1,N2,lambda);
    alpha = alpha + delta(1); d1 = d1 + delta(2); d2 = d2 + delta(3);                              % parameter update

    % make refer data
    y = [Pp(d1,d2) Pp(d1-1,d2) Pp(d1+1,d2) Pp(d1,d2+1) Pp(d1-1,d2+1) Pp(d1+1,d2+1) Pp(d1,d2-1) Pp(d1-1,d2-1) Pp(d1+1,d2-1)].';
    x1 = [d1 d1+1 d1-1 d1 d1+1 d1-1 d1 d1+1 d1-1 ].'; 
    x2 = [d2 d2 d2 d2+1 d2+1 d2+1 d2-1 d2-1 d2-1].';

    cost = cost_function_POC(y,x1,x2,d1,d2,alpha,N1,N2);
    if cost > p_cost
        lambda = 10 * lambda;                                                                                                 % update lambda
    elseif cost > finish_magnitude * init_cost;                                                                         % evaluate cost and not to exit 
        lambda = 0.4 * lambda;                                                                                                % update lambda
        p_cost = cost;                                                                                                               %
    else                                                                                                                                    % finish the roop 
        break;
    end
    

    
    maxroop = maxroop - 1;
end

amp = alpha;
pxx = d1;
pyy = d2;

end

% given [ r , n1 , n2 ] and do least square
function delta = leastsquare_Marquardt(y,x1,x2,d1,d2,alpha,N1,N2,lambda)

lambda = squeeze(lambda);
f = zeros(9,1);
f_a = zeros(9,1);
f_d1 = zeros(9,1);
f_d2 = zeros(9,1);


f = zeros(9,1);
for i = 1:9
    f(i) = alpha/N1/N2 * sin(pi*(x1(i)+d1))*sin(pi*(x2(i)+d2)) / sin(pi*(x1(i)+d1)/N1)*sin(pi*(x2(i)+d2)/N2);
    f_a(i) = f(i)/alpha;
    f_d1(i) = alpha/N1/N2 * sin(pi*(x1(i)+d1)) / sin(pi*(x1(i)+d1)/N1)*sin(pi*(x2(i)+d2)/N2) * pi * ( cos(pi*(x1(i)+d1)) * sin(pi*(x1(i)+d1)/N1) - 1/N1 * sin(pi*(x1(i)+d1)) * cos(pi*(x1(i)+d1)/N1) ) / ( sin(pi*(x1(i)+d1)/N1))^2;
    f_d2(i) = alpha/N1/N2 * sin(pi*(x2(i)+d2)) / sin(pi*(x1(i)+d1)/N1)*sin(pi*(x2(i)+d2)/N2) * pi * ( cos(pi*(x2(i)+d2)) * sin(pi*(x2(i)+d2)/N2) - 1/N2 * sin(pi*(x2(i)+d2)) * cos(pi*(x2(i)+d2)/N2) ) / ( sin(pi*(x2(i)+d2)/N2))^2;
end
E = y -  f;
J =  horzcat(f_a,f_d1,f_d2);

% direct method is not numerically staple :  delta = - ( ( J.' * J + lambda*diag(J.' * J) ) \ J.' ) * E;
% using SVD to avoid this problem
[U S V] = svd(J,0);

sigS = diag(S);
Sigma = sigS./ (sigS.*sigS + lambda^2);
delta = - V*diag(Sigma)*U.' * E;

end

function lambda_0 = init_lambda(y,x1,x2,d1,d2,alpha,N1,N2)

f = zeros(9,1);
f_a = zeros(9,1);
f_d1 = zeros(9,1);
f_d2 = zeros(9,1);


for i = 1:9
    f(i) = alpha/N1/N2 * sin(pi*(x1(i)+d1))*sin(pi*(x2(i)+d2)) / sin(pi*(x1(i)+d1)/N1)*sin(pi*(x2(i)+d2)/N2);
    f_a(i) = f(i)/alpha;
    f_d1(i) = alpha/N1/N2 * sin(pi*(x1(i)+d1)) / sin(pi*(x1(i)+d1)/N1)*sin(pi*(x2(i)+d2)/N2) * pi * ( cos(pi*(x1(i)+d1)) * sin(pi*(x1(i)+d1)/N1) - 1/N1 * sin(pi*(x1(i)+d1)) * cos(pi*(x1(i)+d1)/N1) ) / ( sin(pi*(x1(i)+d1)/N1))^2;
    f_d2(i) = alpha/N1/N2 * sin(pi*(x2(i)+d2)) / sin(pi*(x1(i)+d1)/N1)*sin(pi*(x2(i)+d2)/N2) * pi * ( cos(pi*(x2(i)+d2)) * sin(pi*(x2(i)+d2)/N2) - 1/N2 * sin(pi*(x2(i)+d2)) * cos(pi*(x2(i)+d2)/N2) ) / ( sin(pi*(x2(i)+d2)/N2))^2;
end
E = y -  f;
J =  horzcat(f_a,f_d1,f_d2);

% svd
[U S V] = svd(J,0);
sigS = diag(S);

lambda_0 = max(sigS(:))/100;
end


function cost = cost_function_POC(y,x1,x2,d1,d2,alpha,N1,N2)

f = zeros(9,1);
f_a = zeros(9,1);
f_d1 = zeros(9,1);
f_d2 = zeros(9,1);

for i = 1:9
    f(i) = alpha/N1/N2 * sin(pi*(x1(i)+d1))*sin(pi*(x2(i)+d2)) / sin(pi*(x1(i)+d1)/N1)*sin(pi*(x2(i)+d2)/N2);
    f_a(i) = f(i)/alpha;
    f_d1(i) = alpha/N1/N2 * sin(pi*(x1(i)+d1)) / sin(pi*(x1(i)+d1)/N1)*sin(pi*(x2(i)+d2)/N2) * pi * ( cos(pi*(x1(i)+d1)) * sin(pi*(x1(i)+d1)/N1) - 1/N1 * sin(pi*(x1(i)+d1)) * cos(pi*(x1(i)+d1)/N1) ) / ( sin(pi*(x1(i)+d1)/N1))^2;
    f_d2(i) = alpha/N1/N2 * sin(pi*(x2(i)+d2)) / sin(pi*(x1(i)+d1)/N1)*sin(pi*(x2(i)+d2)/N2) * pi * ( cos(pi*(x2(i)+d2)) * sin(pi*(x2(i)+d2)/N2) - 1/N2 * sin(pi*(x2(i)+d2)) * cos(pi*(x2(i)+d2)/N2) ) / ( sin(pi*(x2(i)+d2)/N2))^2;
end
E = y -  f;
cost = E.' * E;
end