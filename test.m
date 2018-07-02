clear all;
raw=imread('2cm2.png');
raw=sum(double(raw),3);
raw=imresize(raw,[256,256]);
% raw=765-raw;
ctfdeSignificantPix=numel(find(abs(ctfde)>eps(class(ctfde))));

[xsize,ysize]=size(raw);


scatter=raw;
figure;imagesc(scatter);

[xsize,ysize]=size(scatter);
xc=round(xsize/2+1);
yc=round(ysize/2+1);
fourier_temp=fftshift(fft2(scatter));
fourier_temp=fourier_temp.*conj(fourier_temp);
co_scatter=abs(fftshift(ifft2(fourier_temp)));            %自相关
figure;imagesc(co_scatter);

%% subtract the background from the correaltion result
a_mask=zeros(xsize,ysize);                                                %看的是这部分
a_edge=round(xsize*0.15);
a_mask(a_edge:xsize-a_edge+1,a_edge:xsize-a_edge+1)=1;
a=1-a_mask;
a=co_scatter.*a;
a_num=find(a~=0);
a_size=size(a_num(:));
my_mean=sum(sum(a))./max(a_size);
ratio=0.85;
my_mean=my_mean*ratio+max(a(:))*(1-ratio);                      %把边缘的平均光强当作背景值


co_scatter=(co_scatter-my_mean).*(co_scatter-my_mean>0);                  %去背景
a_edge=round(xsize*0.05);
a_mask=zeros(xsize,ysize);
a_mask(a_edge:xsize-a_edge+1,a_edge:xsize-a_edge+1)=1;
a_mask=imgaussfilt(a_mask,3);
% a_maskx=tukeywin(xsize,0.05);
% a_masky=tukeywin(ysize,0.05);
% a_mask=a_maskx*a_masky';
co_scatter=co_scatter.*a_mask;

im_temp=co_scatter;
figure;imagesc(im_temp);pause(0.001);
[xsize,ysize]=size(im_temp);
f_raw=fftshift(fft2(im_temp));                              %S
f_abs=sqrt(abs(f_raw));                                     %S                 
initial_guess=rand([xsize,ysize]);                          %g
a_mask=ones(xsize,ysize).*(co_scatter>=0.01*max(co_scatter(:)));
initial_guess=initial_guess.*a_mask;

phase_temp=exp(1i*angle(fftshift(fft2(initial_guess))));    %G
initial_guess=ifft2(ifftshift(f_abs.*phase_temp));          %g'
initial_guess=abs(initial_guess);                           %初始g


% a_mask=ones(xsize,ysize).*(co_scatter>=0.01*max(co_scatter(:)));
% initial_guess=initial_guess.*a_mask;



im_temp=initial_guess;
b=0.7;                                                      %b的初始值

%% Phase retrieval


%     for jj=1:500                                    
%         gk=fftshift(fft2(im_temp));                                    %G
%         angle_f=angle(gk);                                             %theta
%         g_k=f_abs.*exp(1i.*angle_f);                                   %G'
% %         im_temp_k=abs(ifft2(ifftshift(g_k))).*exp(1i.*imangle);
%         im_temp_k=ifft2(ifftshift(g_k));                               %g'
%     %     mask=(abs(imag(im_temp_k))~=0)+(real(im_temp_k)<0);
%     %     mask=(mask>0.5);
%         mask=(real(im_temp_k)>0).*(abs(imag(im_temp_k))==0);                                      %域的物理筛选条件
%         im_temp=im_temp_k.*mask;                                       %g(k+1） the Error-reduction algorithm.
%          figure(11);imagesc(abs(im_temp));pause(0.001);
%     end
for ii=1:100
    
    
    im_temp_o=im_temp;                                                 %初始g
    
    for jj=1:25
        gk=fftshift(fft2(im_temp));                                    %G
        angle_f=angle(gk);                                             %theta
        g_k=f_abs.*exp(1i.*angle_f);                                   %G'
%         im_temp_k=abs(ifft2(ifftshift(g_k))).*exp(1i.*imangle);
        im_temp_k=ifft2(ifftshift(g_k));                               %g'
    %     mask=(abs(imag(im_temp_k))~=0)+(real(im_temp_k)<0);
    %     mask=(mask>0.5);
        mask=(real(im_temp_k)<0)+(imag(im_temp_k)~=0);                                      %域的物理筛选条件  
        mask=(mask>0.5);
        im_temp=abs(im_temp_k-b.*im_temp.*mask);                          %g(k+1） the HIO algorithm   b从0.7开始以0.95倍衰减若干次
%          figure(12);imagesc(abs(im_temp));pause(0.001);
    end

%  b=b-0.04;       
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
    for jj=1:40                                    
        gk=fftshift(fft2(im_temp));                                    %G
        angle_f=angle(gk);                                             %theta
        g_k=f_abs.*exp(1i.*angle_f);                                   %G'
%         im_temp_k=abs(ifft2(ifftshift(g_k))).*exp(1i.*imangle);
        im_temp_k=ifft2(ifftshift(g_k));                               %g'
    %     mask=(abs(imag(im_temp_k))~=0)+(real(im_temp_k)<0);
    %     mask=(mask>0.5);
        mask=(real(im_temp_k)>0);                                      %域的物理筛选条件
        im_temp=im_temp_k.*mask;                                       %g(k+1） the Error-reduction algorithm.
         figure(11);imagesc(abs(im_temp));pause(0.001);
    end
%% Display results
im_result=abs(im_temp);
[xx,yy]=find(im_result==max(im_result(:)));
im_result=circshift(im_result,[xc-xx(1),yc-yy(1)]);         %移到中心
figure;imagesc(co_scatter);title('correaltion result');
figure;imagesc(im_result);title('reconstructed image');     %输出gk的模大小
% figure;imagesc(co_image);