
load occultation.fla

t=occultation(:,1);
lum=occultation(:,2);

figure(1)
hold off
%plot(t,lum)
%hold on
plot(t,lum,'.r')
xlabel('time, day')
ylabel('Relative luminosity')
title('Simulated transit lightcurve')
print -deps luminosity.eps

figure(2)
nn=size(lum); nl=nn(1); dt=t(2)-t(1)
der=(lum(2:nl)-lum(1:nl-1))./dt;
hold off
%plot(t(1:nl-1),der)
%hold on
plot(t(1:nl-1),der,'.r')
xlabel('time, day')
ylabel('Time derivative of relative luminosity, 1/day')
title('Simulated transit lightcurve derivative')
print -deps derlumin.eps

pause

% store data for future comparison
tp=t;lump=lum;derp=der;nlp=nl;dtp=dt;
clear occultation t lum der
save occsim.mat