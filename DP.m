
%% Intialization
clc
clear all
N=1000;
Alpha=2;
D=2;
Mu0=zeros(1,D);  % Prior mean for cluster Mu
Sig0=1*eye(D);   % Prior covariance matrix for cluster Mu
Sig=0.010*eye(D);    % Prior covariance matrix for data point
Bet=[];
Pi=[];
Mu={};
X=[];
k=0;

%% Dirichlet Process Mixture Model

for n=1:N
    U=rand(1);
    Prev=k;
    I=1;
    if ((sum(Pi))<U)
    while (sum(Pi))<U
        if k~=0
            B=betarnd(1,Alpha,1);
            Pi(Prev+I)=B*prod(1-Bet);
            Bet(Prev+I)=B;
        else
            B=betarnd(1,Alpha,1);
            Pi(Prev+I)=B;
            Bet(Prev+I)=B;
        end
        k=length(Pi);
        I=I+1;
    end
    Curr=length(Pi);
    
    for j=(Prev+1):Curr
        Mu{j}=mvnrnd(Mu0,Sig0,1);
    end
    Z=mvnrnd(Mu{Curr},Sig,1);
    
    else 
       for j=1:length(Pi)
           if sum(Pi(1:j))>U
               MuP=Mu{j};
               break;
           end
       end
       
        Z=mvnrnd(MuP,Sig,1);
    end
    
    scatter(Z(1),Z(2))
    hold on

end
    
for j=1:length(Pi)
    scatter(Mu{j}(1),Mu{j}(2),'filled','d','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0])
    hold on
end
    
    
    
    
            
            
            