%Surface-related multiple elimination SRME is an algorithm process applied to seismic field data.
% SRME using conventional adaptive subtraction
% applied to shot records in a 2-D medium with DIPPING reflector

clear all
Nt=256; % Number of time samples
deltat=0.004; % Time sample distance
T=deltat*Nt;% sampling time
Nx=81; %Number of shot/reciever
fs=1/deltat;% sampling frequency
deltaf=fs/Nt; % frequency sample distance
fmax=60; %maximum frequency in the signal
v1=1500;% Velocity of first subsurface layer
v2=3000;% Velocity of second subsurface layer
Kmax=fmax/v1;
deltax=1/(2*Kmax);
Zlevel=500;
Zrefl=700;% Depth of reflector at X=0
Alpha=0.0872665; % angle of the reflector
R=(v2-v1)/(v2+v1);% Reflection coefficient between two subsurface layers
%p=square(2*pi/1i*k*r)*exp(-1i*k*r);
 is=41;
 xs=(is-1)*deltax;
   for ir=1:Nx
        xr=(ir-1)*deltax;
        p=0;
        h=xr-xs;
       % for imult=1:6
                deltaw=2*pi*deltaf;
               w=(0:deltaw:deltaw*(Nt-1));
                K=w/v1;
                %Zs=Zrefl*cos(Alpha);
                %r=sqrt(((2.*imult.*((Zs+xs.*sin(Alpha))+(h.*sin(Alpha))/2)).^2)+(h.*cos(imult*Alpha)).^2);
                %Teta=asin(h*cos(Alpha)/r);
                %cosphi=cos(Alpha+Teta);
                %p=p+(R^(imult)).*sqrt(2*pi*1i.*w./(v1.*r)).*cosphi.*exp(-1i.*(w/v1).*r);
                r=sqrt((h^2)+Zlevel^2);
                p1=sqrt(2*pi/1i.*(w./v1).*r).*exp(-1i.*(w./v1).*r);
                g1(ir,:)=p1;
       % end
      %%% bandlimiting %%%%
      %%% fmax/deltaf
        for k=1:Nt
            if k<=round(fmax/deltaf)
                M1(k)=(cos((pi*(k-1)/((Nt/2)/2))-pi/2).^2);
            else
                M1(k)=0;
            end
        end
       
        M1=M1.*p1;
        y1=ifft(M1,Nt);
        y1=real(2.*y1./Nt);
        P1(ir,:)=y1;  % Raw data with multiples
        t=(0:(Nt)-1).*deltat;
    end 

%plot the raw data with multiples for middle shot
in11=P1;
in1=in11';
imagesc(((h-500)),t,in1);colormap(gray);
set(gca,'XAxislocation','top','YAxislocation','left','Ydir','reverse');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





Nt=256; % Number of time samples
deltat=0.004; % Time sample distance
T=deltat*Nt;% sampling time
Nx=81; %Number of shot/reciever
fs=1/deltat;% sampling frequency
deltaf=fs/Nt; % frequency sample distance
fmax=60; %maximum frequency in the signal
v1=1500;% Velocity of first subsurface layer
v2=3000;% Velocity of second subsurface layer
Kmax=fmax/v1;
deltax=1/(2*Kmax);
Zlevel1=100;
Zrefl=700;% Depth of reflector at X=0
Alpha=0.3805; % angle of the reflector
R=(v2-v1)/(v2+v1);% Reflection coefficient between two subsurface layers
%p=square(2*pi/1i*k*r)*exp(-1i*k*r);
for is=1:Nx
    xs=(is-1)*deltax;
   for ir=1:Nx
        xr=(ir-1)*deltax;
        p=0;
        h=xr-xs;
       % for imult=1:6
                deltaw=2*pi*deltaf;
               w=(0:deltaw:deltaw*(Nt-1));
                K=w/v1;
                Zlevel2=Zlevel1+xs*atan(Alpha);
                %Zs=Zrefl*cos(Alpha);
                %r=sqrt(((2.*imult.*((Zs+xs.*sin(Alpha))+(h.*sin(Alpha))/2)).^2)+(h.*cos(imult*Alpha)).^2);
                %Teta=asin(h*cos(Alpha)/r);
                %cosphi=cos(Alpha+Teta);
                %p=p+(R^(imult)).*sqrt(2*pi*1i.*w./(v1.*r)).*cosphi.*exp(-1i.*(w/v1).*r);
                r=sqrt((h^2)+Zlevel2^2);
                p2=sqrt(2*pi/1i.*(w./v1).*r).*exp(-1i.*(w./v1).*r);
                
       % end
      %%% bandlimiting %%%%
      %%% fmax/deltaf
        for k=1:Nt
            if k<=round(fmax/deltaf)
                M1(k)=(cos((pi*(k-1)/((Nt/2)/2))-pi/2).^2);
            else
                M1(k)=0;
            end
        end
       
        M2=M1.*p2;
        y1=ifft(M2,Nt);
        y1=real(2.*y1./Nt);
        P2(is,ir,:)=y1;  % Raw data with multiples
        t=(0:(Nt)-1).*deltat;
    end 
end
%plot the raw data with multiples for middle shot
%in11=squeeze(P2(21,:,:));%plot(t,in11)
%in1=in11';
%imagesc(((h-750)),t,in1);colormap(gray);
%set(gca,'XAxislocation','top','YAxislocation','left','Ydir','reverse');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  Exponential time scaling
  alpha=1.3;
        for itime=1:Nt
            P1(:,itime)=exp(-alpha*itime*deltat).*P1(:,itime);
            P2(:,:,itime)=exp(-alpha*itime*deltat).*P2(:,:,itime);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for is=1:Nx

    P3=squeeze(P2(is,:,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
      for ir=1:Nx 
      Daval3(ir,:)=fft(P1(ir,:));
          D125(ir,:)=fft(P3(ir,:));  
      end
 
% Nx monochromatic common shot record are acheived
% multiply the monochromatic shot record to each other as convolotion in time 
for ifreq=1:Nt
    D126=squeeze(D125(:,ifreq));
    Daval4=squeeze(Daval3(:,ifreq));
     DR(:,ifreq)=((-Daval4).*(D126));
    
end

   
      for ir=1:Nx
      DR3(ir,:)=real(ifft(DR(ir,:)));    
      end
  
     

     DR333=0;
        for ir=1:Nx
            DR333=DR3(ir,:)+DR333;
        end
        DRR(is,:)=DR333;
end
% Removing the time scaling    
     for itime=1:Nt
           DR33(:,itime)=exp(alpha*itime*deltat).*DRR(:,itime);
     end
% plot the predicted surface-related multiples
in11=(DR33);%plot(t,in11)
in1=in11';
imagesc(((h)),t,in1);colormap(gray);
set(gca,'XAxislocation','top','YAxislocation','left','Ydir','reverse');


