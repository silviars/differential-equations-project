function dupr_project_62365

%define parameters
 a = sqrt(2/9);
 L = 11;
 tmax = 13;
 x = linspace(0, L);
 t = linspace(0, tmax);

 %define function phi
 function y = phi(x)
   for i = 1:length(x)
     if x(i) >= 1 && x(i) <= 2
       y(i)=5*(1-log(x(i)^2 - 3*x(i) + 3))^3;
     else
       y(i) = 0;
     end
   end
 end

 %define function psi
 function y = psi(x)
   y=sin((pi*x)/11);
 end

 %define function u(x,t)
function y = u(x,t)
 y = 0;
 for k = 0:54
   Xk=sin(((2*k+1)*pi*x)/(2*L));
   Ak=(2/L)*trapz(x,phi(x).*Xk);
   Bk=(4/((2*k+1)*pi*a))*trapz(x,psi(x).*Xk);
   Tk=Ak*cos((2*k+1)*pi*a*t/(2*L))+Bk*sin((2*k+1)*pi*a*t/(2*L));
   y = y + Tk*Xk;
   end
 end

 %generate animation graphics
 for n = 1:length(t)
   plot(x, u(x,t(n)),'r', 'LineWidth', 2)
   axis([0 L -11 11])
   ylabel('u(x,t)');
   grid on
   M(n)=getframe;
 end
 movie(M,1)

 %draw graphics in a single window
 %first stage
 subplot(3,1,1)
 plot(x, u(x, 0), 'r', 'LineWidth', 2)
 axis([0 L -15 15])
 title('t = 0')
 grid on

%mid stage
subplot(3,1,2)
plot(x, u(x, 7), 'r', 'LineWidth', 2)
axis([0 L -11 11])
title('t = 7')
grid on

%final stage
subplot(3,1,3)
plot(x, u(x, tmax), 'r', 'LineWidth', 2)
axis([0 L -11 11])
title('t = 13')
grid on
end
