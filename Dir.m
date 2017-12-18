[X,Y] = meshgrid(0:0.01:1,0:0.01:1);
a=[1/3,1/3,1/3]*0.001;
Z=zeros(length(X(1,:)),length(X(1,:)));
for i=1:length(X(1,:))
    for j=1:length(X(1,:))
        if X(i,j)+Y(i,j)<=1
        Z(i,j)=drchpdf([X(i,j), Y(i,j)],a);
        end
       
    end
end
%v = [-3 -3];
surf(X,Y,Z)