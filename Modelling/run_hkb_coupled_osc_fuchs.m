stoptime=300; %in seconds

% % parameter setting Fuchs 2013 fig 7.3:
% x1_0=0.5;    x2_0=0.5;
% xd1_0=0;   xd2_0=0;
% gamma1=0.7; gamma2=0.7;%linear damping
% eps1=-1;    eps2=-1;% vd pol damping
% ray1=-1;    ray2=-1;%rayleigh damping
% omega1=1.4; omega2=1.4; % eigenfrequency (in rad/s)
% %omega1=2*pi; omega2=2*pi; % eigenfrequency (in rad/s)
% wsq1=-(omega1)^2;    wsq2=-(omega2)^2; % eigenfrequency squared  (in rad/s)
% mu2=1;  % default is 1; mu2 = -1 when competitive accoridng to Kleso 2009 
% eta1= 1;    eta2=0;
% A1=   -0.2;  A2= -0.2;
% B1=  0.2;    B2=0.2;

% %parameter setting Keslo et al 2009:
% x1_0=5;    x2_0=0;
% xd1_0=0;   xd2_0=0;
% gamma1=12.457; gamma2=12.457;%linear damping
% eps1=-0.641;    eps2=-0.641;% vd pol damping
% ray1=-0.00709;    ray2=-0.00709;%rayleigh damping
% omega1=2*pi; omega2=2*pi; % eigenfrequency (in rad/s)
% wsq1=-(omega1)^2;    wsq2=-(omega2)^2; % eigenfrequency squared  (in rad/s)
% mu2=-1;  % default is 1; mu2 = -1 when competitive accoridng to Kleso 2009 
% eta1= 1;    eta2=1;
% A1=   -0.12;  A2= -0.12; % in order to replicate kelso2009 simulation outcomes, both A and B need to be <0 
% B1=  -0.025;    B2=-0.025;
% %B1=  0.1;    B2=0.1;

% my own parameter setting leading to in+anti attraction:
x1_0=0.5;    x2_0=0.5;
xd1_0=0;   xd2_0=0;
gamma1=4; gamma2=4;%linear damping (determining amplitude), must compensate for omega, see HKB 1985
eps1=-8;    eps2=-8;% vd pol damping
ray1=-0.5;    ray2=-0.5;%rayleigh damping
omega1=2*pi; omega2=2*pi; % eigenfrequency (in rad/s)
% omega1=2*pi; omega2=2*pi; % eigenfrequency (in rad/s)
wsq1=-(omega1)^2;    wsq2=-(omega2)^2; % eigenfrequency squared  (in rad/s)
mu2=-0.9; mu1= 1; % default is 1; mu2 = -1 when competitive accoridng to Kleso 2009 
eta1= 1;    eta2=1;
A1=   -0.5;  A2= -0.5;  % increasing this value to e.g. -0.5 leads to monostable regime (see Fuchs 2013)
B1=  0.5;    B2=0.5;
% A1=   -0.2;  A2= -0.2;
% B1=  0.2;    B2=0.2;

%inp_fr=[0.7 0.75 0.8 0.85 0.9];



RP_0=0:pi/8:2*pi; %varying RP-angle in radians with steps of 1/8 pi

parmat=[cos(-RP_0)*x1_0;sin(-RP_0)*x1_0*omega2]';% accoridngly varying initial state of osc 2


%[0 0.1 0.5 0.8];
% figure
% hold on

mycolors=['b','r','g','k','c','m','y','b','r','g','k','c','m','y','b','r','g','k','c','m','y','b','r','g','k','c','m','y'];


RPfig= figure;
hold on
for i = 1:size(parmat,1)
     x2_0=parmat(i,1);
     xd2_0=parmat(i,2);
    %eta1=parmat(i);
    %eta2=eta1;
    %sim('hkb_coupled_osc_fuchs')
    sim('hkb_coupled_osc_kelso2009')


%    figure
%    hold on
%   plot(tout,yout)%, mycolors(i))
%   plot(tout,pos1.signals.values,'r')%,[mycolors(i) ':'] )
%     title(['eta = ' num2str(eta1)])
%  title(['start x xdot =' num2str(x2_0) '  ' num2str(xd2_0)])

    
%      hilph1=cart2pol(yout(:,1),imag(hilbert(yout(:,1))));
%      hilph2=cart2pol(yout(:,2),imag(hilbert(yout(:,2))));
     pks1=peakfind(yout(:,1),0.1,1);
     pks2=peakfind(yout(:,2),0.1,1);
     [imx1, hilph1]=halfcyclehilbert(yout(:,1),pks1);
     [imx2, hilph2]=halfcyclehilbert(yout(:,2),pks2);
     rp=unwrap(hilph1)-unwrap(hilph2);
     
     figure(RPfig)
     if nanmean(rp(end-1000:end))>1.5*pi
         %plot(tout,rp-2*pi,mycolors(i))
         plot(tout,rp-2*pi);plot(tout,rp);plot(tout,rp-4*pi);
     elseif nanmean(rp(end-1000:end))< -0.5*pi
         %plot(tout,rp+2*pi,mycolors(i))
         plot(tout,rp+2*pi);plot(tout,rp+4*pi);plot(tout,rp);
     else
         %plot(tout,rp,mycolors(i))
         plot(tout,rp);plot(tout,rp+2*pi);plot(tout,rp-2*pi);
     end
          figure(3)
          hold on
     if nanmean(rp(end-1000:end))>1.5*pi
         %plot(tout,rp-2*pi,mycolors(i))
         plot(tout,rad2deg(rp-2*pi));plot(tout,rad2deg(rp));plot(tout,rad2deg(rp-4*pi));
     elseif nanmean(rp(end-1000:end))< -0.5*pi
         %plot(tout,rp+2*pi,mycolors(i))
         plot(tout,rad2deg(rp+2*pi));plot(tout,rad2deg(rp+4*pi));plot(tout,rad2deg(rp));
     else
         %plot(tout,rp,mycolors(i))
         plot(tout,rad2deg(rp));plot(tout,rad2deg(rp+2*pi));plot(tout,rad2deg(rp-2*pi));
     end
     
     
    vel1=gradient(yout(:,1),0.001);
    vel2=gradient(yout(:,2),0.001);
  
   

%     figure 
%     plot(yout(:,1),gradient(yout(:,1),0.001)/omega1)
%     hold on
%     plot(yout(:,2),gradient(yout(:,2),0.001)/omega2,'r')
%     plot(x1_0,xd1_0/omega1,'bo')
%     plot(x2_0,xd2_0/omega2,'ro')
%     axis equal

% 
%  figure % to check oscillation amplitudes etc
% plot(yout(:,1));hold on; plot(yout(:,2),'r')
% title(['RP_0='  num2str(RP_0(i)) ' rad'] )




end
