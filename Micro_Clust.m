%% Intialization
clc
clear all
a=50;
q=0.4;
K=nbinrnd(a,q,1); %Number of clusters
while K==0
    K=nbinrnd(a,q,1);
end
F=5; %Number of dimensions
n=0;
i=1;
r=gamrnd(2.2,1,1);
p=betarnd(2,2,1);
j=1;
while(j<=K)
 L=nbinrnd(r,p,1);   
 if L~=0
     n(j)=L;    % Generating N1,N2,...,NK based on the number of clusters
     j=j+1;
 end
    
end 
N=sum(n);           % Total number of points
gam=2*ones(1,F);
for f=1:F
    del(f)=gamrnd(1,2,1);
    for k=1:K
        Theta(f,k)=betarnd(del(f),gam(f),1);
    end
end

S=[];
for j=1:K
    S=[S,ones(1,n(j))*j];
end

j=1;
z=0;
while j<=N
    num=randi([1,length(S)],1);
    z(j)=S(num);
    S(num)=[];
    j=j+1;
end

for w=1:N
    for f=1:F
        x(f,w)=binornd(1,Theta(f,z(w)),1);
    end
end

Iter=1;

%% Algorithm of Micoclustering
% x,f,N,a,q,r,p,Del,Gam are observed
Z=ones(N,1);
%Z(:,1)=1;
R=0;
Iter=1;
while Iter<200
for n0=1:N
    
    P=zeros(1,length(Z(1,:)));
    
    for j=1:length(Z(1,:))
        len=(sum(Z((1:N)~=n0,j)));
        B1=(1+len)*(nbinpdf(len+1,r,p)/nbinpdf(len,r,p));
        B2=prod((x(:,(1:N)~=n0)*Z((1:N)~=n0,j)+del').^x(:,n0).*((1-x(:,(1:N)~=n0))*Z((1:N)~=n0,j)+gam').^(1-x(:,n0)))/(prod((sum(Z((1:N)~=n0,j)))+del+gam));
        P(j)=B1*B2;
    end
   
    
    A1=(length(Z(1,:))+1)*(nbinpdf(length(Z(1,:))+1,a,q)/nbinpdf(length(Z(1,:)),a,q))*(nbinpdf(1,r,p)/nbinpdf(0,r,p));
    A2=prod((del'.^x(:,n0)).*(gam'.^(1-x(:,n0))))/(prod(del+gam));
    
    
    
    P(length(Z(1,:))+1)=A1*A2;
    P=P/sum(P);
    
    H1=mnrnd(1,P);
    H2=1:length(P);
    Index=H2(H1==1);
    if Index<=length(Z(1,:))
        Z(n0,:)=0;
        Z(n0,Index)=1;
    else
        Z(n0,:)=0;
        D=zeros(N,1);
        D(n0)=1;
        Z=[Z,D];
    end
    M=[];
    for j=1:length(Z(1,:))   % removing coulmns of Z which have no point assigned to them
        if(sum(Z(:,j))==0)     
            M=[M,j];
        end
    end
    Z(:,M)=[];
    length(Z(1,:))
    
end
length(Z(1,:))
Iter=Iter+1
R=[R;length(Z(1,:))];
end
%mean(R(2000:5000))  % Number of clusters
        
            
%% Algorithm for DPMM
        
% x,f,N,a,q,r,p,Del,Gam are observed
Z=ones(N,1);
Iter=1;
%Z(:,1)=1;
R=0;
alpha=5;
while Iter<200
for n0=1:N
    
    P=zeros(1,length(Z(1,:)));
    
    for j=1:length(Z(1,:))
        len=(sum(Z((1:N)~=n0,j)));
        B1=(len);
        B2=prod((x(:,(1:N)~=n0)*Z((1:N)~=n0,j)+del').^x(:,n0).*((1-x(:,(1:N)~=n0))*Z((1:N)~=n0,j)+gam').^(1-x(:,n0)))/(prod((sum(Z((1:N)~=n0,j)))+del+gam));
        P(j)=B1*B2;
    end
   
    
    A1=alpha;
    A2=prod((del'.^x(:,n0)).*(gam'.^(1-x(:,n0))))/(prod(del+gam));
    
    
    
    P(length(Z(1,:))+1)=A1*A2;
    P=P/sum(P);
    
    H1=mnrnd(1,P);
    H2=1:length(P);
    Index=H2(H1==1);
    if Index<=length(Z(1,:))
        Z(n0,:)=0;
        Z(n0,Index)=1;
    else
        Z(n0,:)=0;
        D=zeros(N,1);
        D(n0)=1;
        Z=[Z,D];
    end
    M=[];
    for j=1:length(Z(1,:))   % removing coulmns of Z which have no point assigned to them
        if(sum(Z(:,j))==0)     
            M=[M,j];
        end
    end
    Z(:,M)=[];
    length(Z(1,:))
    
end
length(Z(1,:))
Iter=Iter+1
R=[R;length(Z(1,:))];
end
%mean(R(2000:5000))  
    








