clear;
%%%Weibull mortality
maxage=11;
%%%Change Deltat to 1 to obtain the row of values below Deltat=1 in table
%%%1. 
Deltat=0.01;
x=0:Deltat:maxage-1;
c=[.5 1 3];
surv_weib=zeros(length(c),length(x));
b=1;
l_weib=zeros(length(c),length(x));
%survival up to age x
l_weib(3,:)=exp(-1.7*x.^c(1));
l_weib(2,:)=exp(-.5*x.^c(2));
l_weib(1,:)=exp(-.005*x.^c(3));

surv(1,:)=l_weib(1,2:end)./l_weib(1,1:end-1);
surv(2,:)=l_weib(2,2:end)./l_weib(2,1:end-1);
surv(3,:)=l_weib(3,2:end)./l_weib(3,1:end-1);
H_old=zeros(3,1);
H_new=zeros(3,1);
e0=zeros(length(x),1);
e0(1)=1;
for i=1:3
U=diag(surv(i,:),-1);
U(length(x),length(x))=surv(i,end);
temp=ones(1,length(x))-sum(U);
M=diag(temp);
N=(eye(length(x))-U)\eye(length(x));     
eta1=ones(1,length(x))*N;
B=M*N;
etadagger=eta1*B;
eta0=eta1*e0;
H_new(i)=etadagger(1)/eta0
H_old(i)= -l_weib(i,:)*log(l_weib(i,:))'./sum(l_weib(i,:))
end


semilogy(x,l_weib')
legend('.25','1','2')
hold on;
