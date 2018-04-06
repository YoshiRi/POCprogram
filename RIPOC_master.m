% 2016/9/12 Yoshi Ri @ Univ Tokyo
% 2017/8/11 Updated
% RIPOC program 
% input : 2 images
% output : translation , rotation , scaling


%% 画像入力
%AI = rgb2gray(imread('luna1_1.png'));

AI = imread('light.bmp');
BI = imread('dark.bmp');
% 
% AI = imread('normal.bmp');
 BI = imread('light.bmp');
%  BI = AI;

 
%% サイズ決定
[height, width ] = size(AI);
 cy = height/2;
 cx = width/2;

 % Translation, rotation and scaling
 BI = imtranslate(BI,[30, -10]);
 BI = ImageRotateScale(BI,-100,1.2,height,width);

imwrite(uint8(BI),'rotated_dark.png')
%% 窓関数の準備 （画像端の影響を避けるため）
% hannig window and root of hanning window
han_win = zeros(width);
Rhan_win = zeros(width);

% make window 
for i = 1 :height
    for j = 1:width
            han_win(i,j) = 0.25 * (1.0 + cos(pi*abs(cy- i) / height))*(1.0 + cos(pi*abs(cx - j) / width));
            % Root han_win
            Rhan_win(i,j)=abs(cos(pi*abs(cy - i) / height)*cos(pi*abs(cx - j) / width));
            if i >height/8  &&  i<height*7/8
                fi = 1;
            elseif i <= height/8
                fi = (i-1)/height * 8;
            elseif i >= height*7/8
                fi = (height - i)/height *8;
            end
            if j >width/8  &&  j<width*7/8
                fj = 1;
            elseif j <= width/8
                fj = (j-1)/width * 8;
            elseif j >= width*7/8
                fj = (width - j)/width *8;
            end
            trapez_win(i,j)=fi*fj;
    end
end



%% 窓関数（フィルタ）を掛ける(convolute window) 
% IA = Rhan_win .* double(rgb2gray(AI));
 IA = han_win .* double((AI));
 IB = han_win .* double(BI);

 IA = Rhan_win .* double((AI));
 IB = Rhan_win .* double(BI);

% IA = trapez_win .* double((AI));
% IB = trapez_win .* double(BI);

% IA = double((AI));
% IB = double(BI);
% 
%%切り出し
% IA = imcrop(AI,[ cx-width/2,cy-height/2,width-1,height-1]);
% IB = imcrop(BI,[ cx-width/2,cy-height/2,width-1,height-1]); 


%% 2DFFT
A=fft2(IA);
B=fft2(IB);

At = A./abs(A);
Bt = (conj(B))./abs(B);

%% 振幅成分の抽出　＆　対数化
% As=cut_win .* fftshift(log(abs(A)+1));
% Bs=cut_win .* fftshift(log(abs(B)+1));
As= fftshift(log(abs(A)+1));
Bs= fftshift(log(abs(B)+1));

%% Log-Poler Transformation
% need bilinear interpolation

lpcA = zeros(height,width);
lpcB = zeros(height,width);
cx = width / 2;
cy = height / 2;

% cut off val of LPF 
LPmin = width*(1-log2(2*pi)/log2(width));
LPmin = width*(1-log2(4*pi)/log2(width));

%start logplolar 
for i= 0:width-1
        r =power(width,(i)/width);
    for j= 0:height-1
        x=r*cos(2*pi*j/height)+cx;
        y=r*sin(2*pi*j/height)+cy;
        if r < cx * 3 / 5 % in the circle
             x0 = floor(x);
             y0 = floor(y);
             x1 = x0 + 1.0;
             y1 = y0 + 1.0;
            w0=x1-x;
            w1=x-x0;
            h0=y1-y;
            h1=y-y0;
            %　Bilinear補完
            val=As(y0+1,x0+1)*w0*h0 + As(y0+1,x1+1)*w1*h0+ As(y1+1,x0+1)*w0*h1 + As(y1+1,x1+1)*w1*h1;
            %　ほぼ補間でできている低域をCutOffする
            if i > LPmin 
                 lpcA(j+1,i+1)=val;
            else
                 lpcA(j+1,i+1)=0;
            end
            val=Bs(y0+1,x0+1)*w0*h0 + Bs(y0+1,x1+1)*w1*h0+ Bs(y1+1,x0+1)*w0*h1 + Bs(y1+1,x1+1)*w1*h1;
            if i > LPmin 
                 lpcB(j+1,i+1)=val;
            else
                 lpcB(j+1,i+1)=0;
            end
        end
    end
end

%%%% end LogPoler %%%


PA = fft2(lpcA);
PB = fft2(lpcB);
Ap = PA./abs(PA);
Bp = (conj(PB))./abs(PB);
Pp = fftshift(ifft2(Ap.*Bp));


[mm,x]=max(Pp);
[mx,y]=max(mm);
px=y;
py=x(y);
% mado
mg =1;
%% Bilinear補間
% sum_ = Pp(py-1,px-1)+Pp(py,px-1)+Pp(py+1,px-1)+Pp(py-1,px)+Pp(py,px)+Pp(py+1,px)+Pp(py-1,px+1)+Pp(py,px+1)+Pp(py+1,px+1);
% 
% pxx = ( Pp(py-1,px-1)+Pp(py,px-1)+Pp(py+1,px-1) ) * (px-1) + ( Pp(py-1,px)+Pp(py,px)+Pp(py+1,px) ) * px + ( Pp(py-1,px+1)+Pp(py,px+1)+Pp(py+1,px+1) )* (px+1);
% pxx = pxx/sum_;
% 
% pyy = ( Pp(py-1,px-1)+Pp(py-1,px)+Pp(py-1,px+1) ) * (py-1) + ( Pp(py,px-1)+Pp(py,px)+Pp(py,px+1) ) * (py) + ( Pp(py+1,px-1)+Pp(py+1,px)+Pp(py+1,px+1) ) * (py+1);
% pyy= pyy/sum_;
Rect = Pp(py-mg:py+mg,px-mg:px+mg);
vert = sum(Rect,1); horz = sum(Rect,2);
sum11=sum(vert);

pyy=[py-mg:py+mg] * horz /sum11;
pxx=[px-mg:px+mg] * vert' /sum11;

dx = width/2 - pxx + 1;
dy = height/2 - pyy + 1;


%% 回転量には2つのピークが出現する
theta1 = 360 * dy / height;
theta2 = theta1 + 180;
scale = 1/power(width,dx/width)
figure;
mesh(Pp);
ylabel('rotation axis')
xlabel('scaling axis')
zlabel('correlation value')

%% 回転・拡大縮小量 を補正
% 面倒だが角度には2つのパターンがある…
IB_recover1 = ImageRotateScale(IB, theta1,scale,width,height);
IB_recover2 = ImageRotateScale(IB, theta2,scale,width,height);

%% 平行移動量検出 ＆ 回転量決定

IB_R1=fft2(IB_recover1);
IB_R2=fft2(IB_recover2);
IB1p = (conj(IB_R1))./abs(IB_R1);
IB2p = (conj(IB_R2))./abs(IB_R2);

App = A./abs(A);
Pp1 = fftshift(ifft2(App.*IB1p));
Pp2 = fftshift(ifft2(App.*IB2p));

[mm1,x1]=max(Pp1);
[mx1,y1]=max(mm1);
px1=y1;
py1=x1(y1);

[mm2,x2]=max(Pp2);
[mx2,y2]=max(mm2);
px2=y2;
py2=x2(y2);

%% 2種類の回転量についてPOCを行い，ピークが出る方，値が大きいほうを真値とする
if mx1 > mx2
theta = theta1
% bilinear
% sum1 = Pp1(py1-1,px1-1)+Pp1(py1,px1-1)+Pp1(py1+1,px1-1)+Pp1(py1-1,px1)+Pp1(py1,px1)+Pp1(py1+1,px1)+Pp1(py1-1,px1+1)+Pp1(py1,px1+1)+Pp1(py1+1,px1+1);
% 
% pxx1 = ( Pp1(py1-1,px1-1)+Pp1(py1,px1-1)+Pp1(py1+1,px1-1) ) * (px1-1) + ( Pp1(py1-1,px1)+Pp1(py1,px1)+Pp1(py1+1,px1) ) * px1 + ( Pp1(py1-1,px1+1)+Pp1(py1,px1+1)+Pp1(py1+1,px1+1) )* (px1+1);
% pxx1 = pxx1/sum1;
% pyy1 = ( Pp1(py1-1,px1-1)+Pp1(py1-1,px1)+Pp1(py1-1,px1+1) ) * (py1-1) + ( Pp1(py1,px1-1)+Pp1(py1,px1)+Pp1(py1,px1+1) ) * (py1) + ( Pp1(py1+1,px1-1)+Pp1(py1+1,px1)+Pp1(py1+1,px1+1) ) * (py1+1);
% pyy1= pyy1/sum1;

Rect = Pp1(py1-mg:py1+mg,px1-mg:px1+mg);
vert = sum(Rect,1); horz = sum(Rect,2);
sum11=sum(vert);

pyy1=[py1-mg:py1+mg] * horz /sum11;
pxx1=[px1-mg:px1+mg] * vert' /sum11;


% get translation from center
dx = width/2 - pxx1 + 1
dy = height/2 - pyy1 + 1

% show result
f1 = figure;
result = imtranslate(IB_recover1,[-dx, -dy]);
imshow(abs(double(IA)-result),[0 255]);
SaveFigPDF(f1,'sabun_bibun');
f2 = figure;
imshow(result,[0 255]);
SaveFigPDF(f2,'compared_moved');
f3 = figure;
imshow(IA,[0 255]);
SaveFigPDF(f3,'ref_bibun');
f4 = figure;
imshow(IB,[0 255]);
SaveFigPDF(f4,'compared');
figure;
mesh(Pp1);


else
theta = theta2
% sum2 = Pp2(py2-1,px2-1)+Pp2(py2,px2-1)+Pp2(py2+1,px2-1)+Pp2(py2-1,px2)+Pp2(py2,px2)+Pp2(py2+1,px2)+Pp2(py2-1,px2+1)+Pp2(py2,px2+1)+Pp2(py2+1,px2+1);
% 
% pxx2 = ( Pp2(py2-1,px2-1)+Pp2(py2,px2-1)+Pp2(py2+1,px2-1) ) * (px2-1) + ( Pp2(py2-1,px2)+Pp2(py2,px2)+Pp2(py2+1,px2) ) * px2 + ( Pp2(py2-1,px2+1)+Pp2(py2,px2+1)+Pp2(py2+1,px2+1) )* (px2+1);
% pxx2 = pxx2/sum2;
% 
% pyy2 = ( Pp2(py2-1,px2-1)+Pp2(py2-1,px2)+Pp2(py2-1,px2+1) ) * (py2-1) + ( Pp2(py2,px2-1)+Pp2(py2,px2)+Pp2(py2,px2+1) ) * (py2) + ( Pp2(py2+1,px2-1)+Pp2(py2+1,px2)+Pp2(py2+1,px2+1) ) * (py2+1);
% pyy2= pyy2/sum2;

Rect = Pp2(py2-mg:py2+mg,px2-mg:px2+mg);
vert = sum(Rect,1); horz = sum(Rect,2);
sum11=sum(vert);

pyy2=[py2-mg:py2+mg] * horz /sum11;
pxx2=[px2-mg:px2+mg] * vert' /sum11;


dx = width/2 - pxx2 + 1
dy = height/2 - pyy2 + 1

% show result 
f1 = figure;
result = imtranslate(IB_recover2,[-dx, -dy]);
imshow(abs(double(IA)-result),[0 255]);
SaveFigPDF(f1,'sabun_bibun');
f2 = figure;
imshow(result,[0 255]);
SaveFigPDF(f2,'compared_moved');
f3 = figure;
imshow(IA,[0 255]);
SaveFigPDF(f3,'ref_bibun');
f4 = figure;
imshow(IB,[0 255]);
SaveFigPDF(f4,'compared');
figure;
mesh(Pp2);

end

%% evaluation
ref = double(AI);
comp = double(BI);

comp_recover = ImageRotateScale(comp, theta,scale,width,height);
res = imtranslate(comp_recover,[-dx, -dy]);
Error = abs(ref - res);

SSD =  Error(:).' * Error(:)

%% figure
figure();
stiched = [IA,IB;result,abs(double(IA)-result)];
imshow(stiched,[0 255])

stiched = [IA;IB];
imshow(stiched,[0 255])
