function [R,M]=tovf(Dens_1,dens,pres)
c=2.99792458e8; %speed of light
g=6.67E-11; % G constant
i=0;dr=10;
Dens_r(1)=Dens_1;
p0=interp1(dens,pres,Dens_1,'linear','extrap');
p(1)=p0; % center pressure
r(1)=1; %r(1)表示中心小球体半径
m(1)=4*pi/3*Dens_r(1)*r(1)^3; %center point mass,
p(1)=p(1)-2*pi/3*(Dens_r(1)*c^2+p(1))*(Dens_r(1)*c^2+3*p(1))*g/c^4*r(1)^2;

  while i>(-1) %对给定中心密度积分得到M和R
      i=i+1;
      r(i+1)=r(i)+dr;
      k1=-(Dens_r(i)*c^2+p(i))*(m(i)+4*pi*r(i)^3*p(i)/c^2)*g/c^2/(r(i)*(r(i)-2*g/c^2*m(i)));
      l1=4*pi*r(i)^2*Dens_r(i);
      p1=p(i)+dr*k1/2;   %following circle gives the density corresponding the new pressure
      if p1>0 % 如果压强大于0，则
          Dens_r(i)=interp1(pres,dens,p1,'linear','extrap');%  extrap method
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      k2=-(Dens_r(i)*c^2+(p(i)+0.5*dr*k1))*((m(i)+0.5*dr*l1)+4*pi*(r(i)+dr/2)^3*(p(i)+0.5*dr*k1)/c^2)*g/c^2/((r(i)+dr/2)*((r(i)+dr/2)-2*g/c^2*(m(i)+0.5*dr*l1)));
      l2=4*pi*(r(i)+dr/2)^2*Dens_r(i);
      p02=p(i)+dr*k2/2;%following circle gives the density corresponding to new pressure
      if  p02>0
          Dens_r(i)=interp1(pres,dens,p02,'linear','extrap');
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      k3=-(Dens_r(i)*c^2+(p(i)+0.5*dr*k2))*((m(i)+0.5*dr*l2)+4*pi*(r(i)+dr/2)^3*(p(i)+0.5*dr*k2)/c^2)*g/c^2/((r(i)+dr/2)*((r(i)+dr/2)-2*g/c^2*(m(i)+0.5*dr*l2)));
      l3=4*pi*(r(i)+dr/2)^2*Dens_r(i);
      p3=p(i)+dr*k3;%following circle gives the density corresponding the new pressure
      if p3>0
          Dens_r(i)=interp1(pres,dens,p3,'linear','extrap');
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      k4=-(Dens_r(i)*c^2+(p(i)+dr*k3))*((m(i)+dr*l3)+4*pi*(r(i)+dr)^3*(p(i)+dr*k3)/c^2)*g/c^2/((r(i)+dr)*((r(i)+dr)-2*g/c^2*(m(i)+dr*l3)));
      l4=4*pi*(r(i)+dr/2)^2*Dens_r(i);
      p(i+1)=p(i)+1/6*dr*(k1+2*k2+2*k3+k4);
      m(i+1)=m(i)+1/6*dr*(l1+2*l2+2*l3+l4);
      % 4 rank runge-Kutta integral method
      
      if p(i+1)>0
          Dens_r(i+1)=interp1(pres,dens,p(i+1),'linear','extrap');
      else
          j=i+1;
          i=-1;
          Dens_r(j)=Dens_r(j-1);
          p(j)=0; %使星体表面压强为准确零.
      end
      
  end
  M=m(j-1);
  R=r(j);
end

