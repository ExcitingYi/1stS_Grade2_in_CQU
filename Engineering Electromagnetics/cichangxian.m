clc;clear;
x=[0:0.001:1];
z=0.*x;
y = 1.+z;
%画矩形线圈。
plot3(x,y,z,'k');
hold on;
plot3(x,z,z,'k');
plot3(z,x,z,'k')
plot3(y,x,z,'k');
%从下往上，一条一条画，先画左半边。
for j = -0.001:-1:-5.001    %从z = j开始，从下往上画。
    x = 0.001;      %从x=y=0.001开始画。
    y = 0.001;
    t = zeros(1,2*j*1000+1);    %t作为坐标x和y，开始画磁力线
    z = [j:0.001:-j];    %重新定义z
    k = 1;  %计数。
    for i = j:0.001:-j     
        %因为是要求比例，Bz/Bxy = dz/dxy,所以u0,I和4派都省略了。
        B1x = (1/sqrt(x^2+i^2))*(y/sqrt(x^2+y^2+i^2) + (1-y)/sqrt(x^2+(1-y)^2+i^2))*(i/sqrt(x^2+i^2));
        B2x = (1/sqrt((1-x)^2+i^2))*(y/sqrt((1-x)^2+y^2+i^2) + (1-y)/sqrt((1-x)^2+(1-y)^2+i^2))*(-i/sqrt((1-x)^2+i^2));
        B3y = (1/sqrt((1-y)^2+i^2))*(x/sqrt(x^2+(1-y)^2+i^2) + (1-x)/sqrt((1-x)^2+(1-y)^2+i^2))*(-i/sqrt((1-y)^2+i^2));
        B4y = (1/sqrt(y^2+i^2))*(x/sqrt(x^2+y^2+i^2) + (1-x)/sqrt((1-x)^2+y^2+i^2))*(i/sqrt(y^2+i^2));
        a = sqrt(2)/2*(B1x+B2x+B3y+B4y);
        B1z = (1/sqrt(x^2+i^2))*(y/sqrt(x^2+y^2+i^2) + (1-y)/sqrt(x^2+(1-y)^2+i^2))*(-x/sqrt(x^2+i^2));
        B2z = (1/sqrt((1-x)^2+i^2))*(y/sqrt((1-x)^2+y^2+i^2) + (1-y)/sqrt((1-x)^2+(1-y)^2+i^2))*(-(1-x)/sqrt((1-x)^2+i^2));
        B3z = (1/sqrt((1-y)^2+i^2))*(x/sqrt(x^2+(1-y)^2+i^2) + (1-x)/sqrt((1-x)^2+(1-y)^2+i^2))*(-(1-y)/sqrt((1-y)^2+i^2));
        B4z = (1/sqrt(y^2+i^2))*(x/sqrt(x^2+y^2+i^2) + (1-x)/sqrt((1-x)^2+y^2+i^2))*(-y/sqrt(y^2+i^2));
        b = B1z + B2z + B3z + B4z;
        xy = 0.001 * a / b;
        x = x+xy/sqrt(2);
        y = x;
        t(k) = x;
        k = k+1;
    end
    plot3(t,t,z,'b');
end
%画右半边。
for j = -0.001:-1:-5.001    %从z = j开始，从下往上画。
    x = 0.999;      %从x=y=0.001开始画。
    y = 0.999;
    t = zeros(1,2*j*1000+1);    %t作为坐标x和y，开始画磁力线。
    z = [j:0.001:-j];    %重新定义z
    k = 1;  %计数。
    for i = j:0.001:-j     
        %因为是要求比例，Bz/Bxy = dz/dxy,所以u0,I和4派都省略了。
        B1x = (1/sqrt(x^2+i^2))*(y/sqrt(x^2+y^2+i^2) + (1-y)/sqrt(x^2+(1-y)^2+i^2))*(i/sqrt(x^2+i^2));
        B2x = (1/sqrt((1-x)^2+i^2))*(y/sqrt((1-x)^2+y^2+i^2) + (1-y)/sqrt((1-x)^2+(1-y)^2+i^2))*(-i/sqrt((1-x)^2+i^2));
        B3y = (1/sqrt((1-y)^2+i^2))*(x/sqrt(x^2+(1-y)^2+i^2) + (1-x)/sqrt((1-x)^2+(1-y)^2+i^2))*(-i/sqrt((1-y)^2+i^2));
        B4y = (1/sqrt(y^2+i^2))*(x/sqrt(x^2+y^2+i^2) + (1-x)/sqrt((1-x)^2+y^2+i^2))*(i/sqrt(y^2+i^2));
        a = sqrt(2)/2*(B1x+B2x+B3y+B4y);
        B1z = (1/sqrt(x^2+i^2))*(y/sqrt(x^2+y^2+i^2) + (1-y)/sqrt(x^2+(1-y)^2+i^2))*(-x/sqrt(x^2+i^2));
        B2z = (1/sqrt((1-x)^2+i^2))*(y/sqrt((1-x)^2+y^2+i^2) + (1-y)/sqrt((1-x)^2+(1-y)^2+i^2))*(-(1-x)/sqrt((1-x)^2+i^2));
        B3z = (1/sqrt((1-y)^2+i^2))*(x/sqrt(x^2+(1-y)^2+i^2) + (1-x)/sqrt((1-x)^2+(1-y)^2+i^2))*(-(1-y)/sqrt((1-y)^2+i^2));
        B4z = (1/sqrt(y^2+i^2))*(x/sqrt(x^2+y^2+i^2) + (1-x)/sqrt((1-x)^2+y^2+i^2))*(-y/sqrt(y^2+i^2));
        b = B1z + B2z + B3z + B4z;
        xy = 0.001 * a / b;
        x = x+xy/sqrt(2);
        y = x;
        t(k) = x;
        k = k+1;
    end
    plot3(t,t,z,'b');
end

view(30,45);
% function [a,b]=B(x,y,z)
%     %因为是要求比例，Bz/Bxy = dz/dxy,所以u0,I和4派都省略了。
%     B1x = (1/sqrt(x^2+z^2))*(y/sqrt(x^2+y^2+z^2) + (1-y)/sqrt(x^2+(1-y)^2+z^2))*(z/sqrt(x^2+z^2));
%     B2x = (1/sqrt((1-x)^2+z^2))*(y/sqrt((1-x)^2+y^2+z^2) + (1-y)/sqrt((1-x)^2+(1-y)^2+z^2))*(-z/sqrt((1-x)^2+z^2));
%     B3y = (1/sqrt((1-y)^2+z^2))*(x/sqrt(x^2+(1-y)^2+z^2) + (1-x)/sqrt((1-x)^2+(1-y)^2+z^2))*(-z/sqrt((1-y)^2+z^2));
%     B4y = (1/sqrt(y^2+z^2))*(x/sqrt(x^2+y^2+z^2) + (1-x)/sqrt((1-x)^2+y^2+z^2))*(z/sqrt(y^2+z^2));
%     B1z = (1/sqrt(x^2+z^2))*(y/sqrt(x^2+y^2+z^2) + (1-y)/sqrt(x^2+(1-y)^2+z^2))*(-x/sqrt(x^2+z^2));
%     B2z = (1/sqrt((1-x)^2+z^2))*(y/sqrt((1-x)^2+y^2+z^2) + (1-y)/sqrt((1-x)^2+(1-y)^2+z^2))*(-(1-x)/sqrt((1-x)^2+z^2));
%     B3z = (1/sqrt((1-y)^2+z^2))*(x/sqrt(x^2+(1-y)^2+z^2) + (1-x)/sqrt((1-x)^2+(1-y)^2+z^2))*(-(1-y)/sqrt((1-y)^2+z^2));
%     B4z = (1/sqrt(y^2+z^2))*(x/sqrt(x^2+y^2+z^2) + (1-x)/sqrt((1-x)^2+y^2+z^2))*(-y/sqrt(y^2+z^2));
%     a = sqrt(2)/2*(B1x+B2x+B3y+B4y);
%     b = B1z + B2z + B3z + B4z;
% end