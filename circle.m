function [ A ] = circle( r,xsize,ysize )
%r:radius of curvature;
%xsize&ysize: to normalize the circle
%
A=zeros(xsize,ysize);
for ii=1:xsize
    for jj=1:ysize
        if ((ii-xsize/2-1)/xsize*r)^2+((jj-ysize/2-1)/ysize*r)^2<=1
            A(ii,jj)=1;%sin(sqrt((ii-xsize/2-1)^2+(jj-ysize/2-1)^2)/sqrt((xsize/r)^2+(ysize/r)^2)*pi)/(sqrt((ii-xsize/2-1)^2+(jj-ysize/2-1)^2)/sqrt((xsize/r)^2+(ysize/r)^2)*pi);
        end
    end
end
A(floor(xsize/2+1),floor(ysize/2+1))=1;

end

