raw=imread('huan1.png');

raw=sum(double(raw),3);

% raw=imresize(raw(127-64:127+64,127-64:127+64),[256,256]);
raw=imresize(raw,[256,256]);
% raw_min=min(min(raw));
% raw=(raw-200).*(raw-200>0);
% raw_temp=imresize(raw,[512,512]);
% raw=zeros(1024,1024);
% raw(256:255+512,256:255+512)=raw_temp;
[xsize,ysize]=size(raw);

[Y,X]=meshgrid(1:ysize,1:xsize);
xc=round(xsize/2+1);
yc=round(ysize/2+1);
yr=Y-yc;
xr=X-xc;

lamda=670;
pixelsize=(5500/1)/50;      %这啥啊      
NA=1.49;
% NA=0.25;

%% Generation of the PSF with Besselj.
R=sqrt((xr).^2+(yr).^2);
pixelnum=xsize;
rpixel=NA*pixelnum*pixelsize/lamda;          %与下面圆的半径有关，但这表达式什么意思
[M1,N1]=meshgrid(1:pixelnum,1:pixelnum);
ctfde=ones(pixelnum,pixelnum).*(((N1-(pixelnum+1)/2)/rpixel).^2+((M1-(pixelnum+1)/2)/rpixel).^2<=1);    %画了个圆
ctfdeSignificantPix=numel(find(abs(ctfde)>eps(class(ctfde))));
ifftscalede=numel(ctfde)/ctfdeSignificantPix;       %矩阵中圆的比例的倒数

apsfde=fftshift(ifft2(ifftshift(ctfde)));
ipsfde=ifftscalede*abs(apsfde).^2;
OTFde=fftshift(fft2(ifftshift(ipsfde)));            %？处理之后的圆（低通滤波器）

figure;imagesc(abs(OTFde));title('Generation of PSF');
cutoff=OTFedgeF(abs(OTFde));
cutoff=cutoff+round(cutoff*0.05);                  %cutoff的意义是什么
f_raw=fftshift(fft2(raw));
raw=ifft2(ifftshift(f_raw.*OTFde));                %原图像进行了低通滤波

fmask=circle(1.1*(xsize+ysize)./(2*cutoff),xsize,ysize);
fmask=imgaussfilt(fmask,10);                       %又花了个圆，半径由cutoff决定
                                                   %圆是模拟成像透镜的低通滤波么，为什么要两个圆，而且还进行了模糊

my_angle=pi/3;
kx=cos(my_angle);
ky=sin(my_angle);
r=6;
phase=(xr.*kx+yr.*ky).*2*pi/r;
pattern_temp=0.5-0.5*cos(phase);
scatter = random('Normal',0,2*pi,xsize,ysize);
% scatter=abs(scatter);
% s_mean=mean2(scatter);
% s_mean=s_mean.*0.25;
% scatter=scatter-s_mean;
pmask=ctfde.*exp(1i*scatter);
apsf=fftshift(ifft2(ifftshift(pmask)));
ipsf=ifftscalede*abs(apsf).^2;
ipsf=ipsf-min(ipsf(:));
OTF_scatter=fftshift(fft2(ifftshift(ipsf)));                     %随机psf的生成
% pmask=pmask./abs(pmask);
i_f=fftshift(fft2(pattern_temp));
% test=i_f.*pmask;
% im=ifft2(ifftshift(test));
% figure;imagesc(im);
% figure;imagesc(pattern_temp);
im_test=fftshift(fft2(raw)).*OTF_scatter;
% im_test=abs(ifft2(ifftshift(im_test.*OTFo)));
im_test=abs(ifft2(ifftshift(im_test)));
im_test=im_test-min(min(im_test));
im_test=im_test./max(im_test(:))*xsize*ysize;
% raw=im_test;
% figure;imagesc(im_test);

% figure;imagesc(raw);

% co_scatter=zeros(xsize,ysize);
% co_image=zeros(xsize,ysize);
% % co_raw=zeros(xsize,ysize);
% % scatter=scatter-mean2(scatter);
% for ii=1:xsize
%     x_t=mod(ii+xc-2,xsize);
%     x_t=x_t+1;
%     for jj=1:ysize
%         y_t=mod(jj+yc-2,ysize);
%         y_t=y_t+1;
% %         co_scatter(x_t,y_t)=sum(sum(scatter.*circshift(scatter,[ii-1,jj-1])));
%         co_image(x_t,y_t)=sum(sum(im_test.*circshift(im_test,[ii-1,jj-1])));
% %         co_raw(x_t,y_t)=sum(sum(raw.*circshift(raw,[ii-1,jj-1])));
%     end
% end
% scatter=imresize(im_test,2,'bicubic');
scatter=im_test;

% scatter=raw;

[xsize,ysize]=size(scatter);
fourier_temp=fftshift(fft2(scatter)).*fmask;
fourier_temp=fourier_temp.*conj(fourier_temp);
% co_scatter=ifft2(ifftshift(fourier_temp));
co_scatter=abs(fftshift(ifft2(fourier_temp)));            %自相关，我把这里的fft改成了ifft
figure;imagesc(co_scatter);pause(0.001); title('Raw autocorrelation'); %%% Autocorrelation


%% Subtracts background from the correlation result

%%% We will calculate the mean in two different ways and take the max value
my_mean = zeros(1,2);
edge_temp = 8;

%%% 1) Gets mean value of corners
%%% Gets the four corners of co_scatter
scatter_corners = co_scatter([1:edge_temp, end-edge_temp+1: end], ...
                             [1:edge_temp, end-edge_temp+1: end]);
my_mean(1) = mean2(scatter_corners);

%%% ---Original code---
% a=co_scatter(1:edge_temp,1:edge_temp);
% a_mean=mean2(a);
% a=co_scatter(xsize-edge_temp+1:xsize,1:edge_temp);
% a_mean=a_mean+mean2(a);
% a=co_scatter(1:edge_temp,ysize-edge_temp+1:ysize);
% a_mean=a_mean+mean2(a);
% a=co_scatter(xsize-edge_temp+1:xsize:xsize,ysize-edge_temp+1:ysize);
% a_mean=a_mean+mean2(a);
% my_mean=a_mean/4;

%%% 2) Gets mean value of a frame and weights with max value element in it
weight = 0.85;

edge_temp = round(xsize*0.08);
mask_temp = zeros(xsize, ysize);
mask_temp(edge_temp:xsize-edge_temp+1,edge_temp:xsize-edge_temp+1) = 1;
ones_frame = 1-mask_temp;
value_frame = co_scatter.*ones_frame;
frame_count = numel(find(value_frame));
my_mean(2) = weight*sum(sum((value_frame))/frame_count) + ...
             (1-weight)*max(value_frame(:));

%%% ---Original code---
% a_mask=zeros(xsize,ysize);                                                
% edge_temp=round(xsize*0.08);
% a_mask(edge_temp:xsize-edge_temp+1,edge_temp:xsize-edge_temp+1)=1;
% a=1-a_mask;
% a=co_scatter.*a;
% a_num=find(a~=0);
% a_size=size(a_num(:));
% my_mean(2)=sum(sum(a))./max(a_size);
% ratio=0.85;
% my_mean(2)=my_mean(2)*ratio+max(a(:))*(1-ratio);    %%% Takes most intensive pixel into account
%my_mean=max(my_mean); 

max_mean = max(my_mean);


co_scatter=(co_scatter-max_mean).*(co_scatter-max_mean>0);                 
edge_temp=round(xsize*0.15);
mask_temp=zeros(xsize,ysize);
mask_temp(edge_temp:xsize-edge_temp+1,edge_temp:xsize-edge_temp+1)=1;
gauss_filter=imgaussfilt(mask_temp,3); %%% SD = 3
co_scatter=co_scatter.*gauss_filter;% force edges to 0                          


% PSFd = real(fftshift( ifft2(fftshift(abs(OTFde).^3)) ));
% PSFd = PSFd/max(max(PSFd));
% PSFd = PSFd/sum(sum(PSFd));
% h = 30;
% PSFe = PSFd(xc-h+1:xc+h,yc-h+1:yc+h);
% 
% co_scatter=edgetaper(co_scatter,PSFe);


im_temp=co_scatter;
figure;imagesc(im_temp);pause(0.001); title('Correlation result 2');

return; %%% Stops program (temp)

%% Phase retrieval

% im_temp=(co_scatter./max(max(co_scatter)))*10^4;
% im_temp=(im_temp-50).*(im_temp-50>0);
[xsize,ysize]=size(im_temp);
f_raw=fftshift(fft2(im_temp));                              %R
f_abs=sqrt(abs(f_raw));                                     %S                 
initial_guess=rand([xsize,ysize]);                          %g
phase_temp=exp(1i*angle(fftshift(fft2(initial_guess))));    %G
initial_guess=ifft2(ifftshift(f_abs.*phase_temp));          %g'
initial_guess=abs(initial_guess);                           %?
% my_max=max(max(initial_guess));       
% initial_guess(xc-5:xc+5,yc-5:yc+5)=my_max;

% fmask=circle((xsize+ysize)./(2*cutoff),xsize,ysize);
% f_abs=f_abs.*fmask;

f_initial=fftshift(fft2(initial_guess));
f_initial=f_initial.*fmask;                                 %最前面的fmask，低通滤波
im_temp=abs(ifft2(ifftshift(f_initial)));                   %?未知处理后的的g',变成了初始下一部分的初始g
% im_temp=ones(xsize,ysize);
b=0.7;                                                      %b的初始值
% angle_f=0.5*angle(f_raw);
% im_temp=abs(fftshift(fft2(f_temp)));
% angle_f=rand([xsize,ysize]);

for ii=1:100

    
    im_temp_o=im_temp;                                                 %初始g
    
    
    for jj=1:5
        imangle=angle(im_temp);                                        %无效
        gk=fftshift(fft2(im_temp));                                    %G
        angle_f=angle(gk);                                             %theta
        g_k=f_abs.*exp(1i.*angle_f);                                   %G'
%         im_temp_k=abs(ifft2(ifftshift(g_k))).*exp(1i.*imangle);
        im_temp_k=ifft2(ifftshift(g_k));                               %g'
    %     mask=(abs(imag(im_temp_k))~=0)+(real(im_temp_k)<0);
    %     mask=(mask>0.5);
        mask=(real(im_temp_k)>0);                                      %域的物理筛选条件
        im_temp=im_temp_k.*mask;                                       %g(k+1） the Error-reduction algorithm.
%         figure(11);imagesc(abs(im_temp));
    end
    
    for jj=1:25
%%%     imangle=angle(im_temp);                                        %无效
        gk=fftshift(fft2(im_temp));                                    %G
        angle_f=angle(gk);                                             %theta
        g_k=f_abs.*exp(1i.*angle_f);                                   %G'
%         im_temp_k=abs(ifft2(ifftshift(g_k))).*exp(1i.*imangle);
        im_temp_k=ifft2(ifftshift(g_k));                               %g'
    %     mask=(abs(imag(im_temp_k))~=0)+(real(im_temp_k)<0);
    %     mask=(mask>0.5); 
        mask=(real(im_temp_k)<0);                                      %域的物理筛选条件   
        im_temp=(im_temp_k-b.*im_temp.*mask);                          %g(k+1） the HIO algorithm   b从0.7开始以0.95倍衰减若干次
%         figure(11);imagesc(abs(im_temp));
    end

        
    MSE(1,ii)=sum(sum(abs(im_temp_o-im_temp)));                        %目前g（k）与初始g的差的总和
    if ii>1&&MSE(1,ii)>MSE(1,ii-1)                                     %如果比前一个循环大 ？？
        b=b*0.95;                                                      %b衰减
    end
        figure(11);imagesc(abs(im_temp));pause(0.001);
%     temp_im1=im_temp;
%     temp_FT=fftshift(fft2(im_temp));
%     lowFT=OTFo.*temp_FT;
%     lowtemp_im=ifft2(ifftshift(OTFo.*temp_FT));
%     lowtemp_im=imlowseq(:,:,i).*exp(1i.*angle(lowtemp_im));
%     lowtemp_FT=fftshift(fft2(lowtemp_im));
%     temp_FT=temp_FT+conj(fmask)./(max(max((abs(fmask)).^2))).*(lowtemp_FT-lowFT);
%     temp_im=ifft2(ifftshift(temp_FT));
%     un_pattern=un_pattern+him./(max(max((abs(him)).^2))+eps).*(temp_im-temp_im1);  
%     him=him+(un_pattern).*(temp_im-temp_im1)./(max(max(un_pattern))).^2;
%     un_pattern1(p:q,p:q)=un_pattern;     
%     himFT=fftshift(fft2(him));
    
end
title('Phase retrieval');
% for ii=1:100
%     
%     
%     im_temp_o=im_temp;                                                 %初始g
%     
%     for jj=1:25
%         gk=fftshift(fft2(im_temp));                                    %G
%         angle_f=angle(gk);                                             %theta
%         g_k=f_abs.*exp(1i.*angle_f);                                   %G'
% %         im_temp_k=abs(ifft2(ifftshift(g_k))).*exp(1i.*imangle);
%         im_temp_k=ifft2(ifftshift(g_k));                               %g'
%     %     mask=(abs(imag(im_temp_k))~=0)+(real(im_temp_k)<0);
%     %     mask=(mask>0.5);
%         mask=(real(im_temp_k)<0)+(imag(im_temp_k)~=0);                                      %域的物理筛选条件  
%         mask=(mask>0.5);
%         im_temp=abs(im_temp_k-b.*im_temp.*mask);                          %%% FEL? %g(k+1） the HIO algorithm   b从0.7开始以0.95倍衰减若干次
% %          figure(12);imagesc(abs(im_temp));pause(0.001);
%     end
% 
% %  b=b-0.04;       
%     MSE(1,ii)=sum(sum(abs(im_temp_o-im_temp)));                        %目前g（k）与初始g的差的总和
%     if ii>1&&MSE(1,ii)>MSE(1,ii-1)                                     %如果比前一个循环大 ？？
%         b=b*0.95;                                                      %b衰减
%     end
%         figure(11);imagesc(abs(im_temp));pause(0.001);
% %     temp_im1=im_temp;
% %     temp_FT=fftshift(fft2(im_temp));
% %     lowFT=OTFo.*temp_FT;
% %     lowtemp_im=ifft2(ifftshift(OTFo.*temp_FT));
% %     lowtemp_im=imlowseq(:,:,i).*exp(1i.*angle(lowtemp_im));
% %     lowtemp_FT=fftshift(fft2(lowtemp_im));
% %     temp_FT=temp_FT+conj(fmask)./(max(max((abs(fmask)).^2))).*(lowtemp_FT-lowFT);
% %     temp_im=ifft2(ifftshift(temp_FT));
% %     un_pattern=un_pattern+him./(max(max((abs(him)).^2))+eps).*(temp_im-temp_im1);  
% %     him=him+(un_pattern).*(temp_im-temp_im1)./(max(max(un_pattern))).^2;
% %     un_pattern1(p:q,p:q)=un_pattern;     
% %     himFT=fftshift(fft2(him));
%     
% end
%     for jj=1:40                                    
%         gk=fftshift(fft2(im_temp));                                    %G
%         angle_f=angle(gk);                                             %theta
%         g_k=f_abs.*exp(1i.*angle_f);                                   %G'
% %         im_temp_k=abs(ifft2(ifftshift(g_k))).*exp(1i.*imangle);
%         im_temp_k=ifft2(ifftshift(g_k));                               %g'
%     %     mask=(abs(imag(im_temp_k))~=0)+(real(im_temp_k)<0);
%     %     mask=(mask>0.5);
%         mask=(real(im_temp_k)>0);                                      %域的物理筛选条件
%         im_temp=im_temp_k.*mask;                                       %g(k+1） the Error-reduction algorithm.
%          figure(11);imagesc(abs(im_temp));pause(0.001);
%     end
%% Display results
im_result=abs(im_temp);
[xx,yy]=find(im_result==max(im_result(:)));
im_result=circshift(im_result,[xc-xx(1),yc-yy(1)]);         %移到中心
figure;imagesc(co_scatter);title('Correlation result');
figure;imagesc(im_result);title('Reconstructed image');     %输出gk的模大小
% figure;imagesc(co_image);