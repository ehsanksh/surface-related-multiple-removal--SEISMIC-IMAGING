%Surface-related multiple elimination SRME is an algorithm process applied to seismic field data.
% SRME using conventional adaptive subtraction
% applied to shot records in a 2-D medium with reflector

clear all
clc
Nt=1024;
deltat=0.004;
T=deltat*Nt;
Nx=6;
fs=1/deltat;
deltaf=fs/Nt;
fmax=60;
v1=1500;
v2=3000;
Kmax=fmax/v1;
deltax=1/(2*Kmax);
Zrefl=300;
Alpha=0.0872665;
R=(v2-v1)/(v2+v1);
for is=1:Nx
    for ir=1:Nx
        xs=(is-1)*deltax;
        xr=(ir-1)*deltax;
        p=0;
        h=xr-xs;
        for imult=1:6
            for k=1:Nt/2
                deltaw=2*pi*deltaf;
                w=(k-1).*deltaw;
                w=(0:deltaw:deltaw*(Nt/2));
                K=w/v1;
                Zs=Zrefl*cos(Alpha);
                r=sqrt(((2.*imult.*((Zs+xs.*sin(Alpha))+(h.*sin(Alpha))/2)).^2)+(h.*cos(imult*Alpha)).^2);
                Teta=asin(h*cos(Alpha)/r);
                cosphi=cos(Alpha+Teta);
                p=p+(R^(imult)).*sqrt(2*pi*1i.*w./(v1.*r)).*cosphi.*exp(-1i.*(w/v1).*r);
            end
        end
      %%% bandlimiting %%%%
      %%% fmax/deltaf=246
        for k=1:(Nt/2)+1
            if k<=round(fmax/deltaf)
                M1(k)=(cos((pi*(k-1)/((Nt/2)/2))-pi/2).^2);
            else
                M1(k)=0;
            end
        end
       
        M=M1.*p;
        y1=ifft(M);
        y1=ifft(M,Nt);
        y1=y1(1:(Nt/2));
        y1=real(2.*y1./Nt);
        P(is,ir,:)=y1;  % Raw data with multiples
        t=(0:(Nt/2)-1).*deltat;
    end 
end
% plot the raw data with multiples for middle shot
%in11=squeeze(P(41,:,:));%plot(t,in11)
 %in1=in11';
%imagesc(h,t,in1);colormap(gray);
%set(gca,'XAxislocation','top','YAxislocation','left','Ydir','reverse'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%  Exponential time scaling
  alpha=1.3;
        for itime=1:Nt/2
            P(:,:,itime)=exp(-alpha*itime*deltat).*P(:,:,itime);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D12=P;
for q=1:4 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     for is=1:Nx
      for ir=1:Nx 
      Daval3(is,ir,:)=fft(P(is,ir,:));
          D125(is,ir,:)=fft(D12(is,ir,:));  
      end
end 
% Nx chromatic common shot record are acheived
%%% multiply the monochromatic shot record to each other as convolotion 
for ifreq=1:Nt/2
    D126=squeeze(D125(:,:,ifreq));
    Daval4=squeeze(Daval3(:,:,ifreq));
     DR(:,:,ifreq)=((-Daval4)*(D126));
    
end
   for is=1:Nx
      for ir=1:Nx
      DR3(is,ir,:)=real(ifft(DR(is,ir,:)));
          
      end
  end
     % Removing the time scaling    
     for itime=1:Nt/2
           DR33(:,:,itime)=exp(alpha*itime*deltat).*DR3(:,:,itime);
        end   

% plot the predicted surface multiples
%in11=squeeze(DR33(41,:,:));%plot(t,in11)
 %imagesc(h,t,in1);colormap(gray);
%set(gca,'XAxislocation','top','YAxislocation','left','Ydir','reverse'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx=Nx;
lenf=7; 
nta=Nt/2;
ntb=Nt/2;
l=lenf;
nxb=nx;
nxa=nx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for is=1:Nx
in1(:,:)=DR33(is,:,:);ref1(:,:)=P(is,:,:);
in=in1';
ref=ref1';

nrot=floor((lenf-1)/2);
zz=zeros(nrot,nx);

        a=in;
        b=in;
        %[nt,nxa]=size(a);
%[nt,nxb]=size(b);
nt=max(nta,ntb);
aa=zeros(nt,nx);
bb=zeros(nt,nx);
c=zeros(l,nx);
aa(1:nt,:)=a;
bb(1:nt,:)=b;

for j=1:l
    c(j,:)=sum(aa(j:nt,:).*bb(1:nt+1-j,:));  
end
%c=c1';
 xx=sum(c,2);
        %xx=sum(correv(in,in,lenf)')';
        a=[zz;ref];
        b=[in;zz];
        [nt1,nxa1]=size(a);
[nt1,nxb1]=size(b);
aa=zeros(nt1,nxa1);
bb=zeros(nt1,nxb1);
c=zeros(l,nxa1);
aa(1:nt1,:)=a;
bb(1:nt1,:)=b;
lenf=l;
for j=1:l
    c(j,:)=sum(aa(j:nt1,:).*bb(1:nt1+1-j,:));  
end

%cc=c2';
yx=sum(c,2);

    Txx=toeplitz(xx);
    Tmax=max(max(Txx));
    Tstab=eps*eps*Tmax*eye(lenf);
    ff=(Tstab+Txx)\yx;
    out=conv2(in,ff,'same');
   out=out';
   D124=ref1-out;
   D12(is,:,:)=D124;
   out1(is,:,:)=out;
end
   D12=P-out1;
end
   in11=squeeze(D12(3,:,:));%plot(t,in1)
 in1=in11';
imagesc(h,t,in1);colormap(gray);
set(gca,'XAxislocation','top','YAxislocation','left','Ydir','reverse');           
   
            
    
    
        