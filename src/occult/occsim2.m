

load occultation.fla
load occsim.mat

t=occultation(:,1);
lum=occultation(:,2);

figure(1)
hold off
%plot(t,lum)
%hold on
plot(t,lum,'.r') % red is new
hold on
plot(tp,lump,'.k') % black is old
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
plot(t(1:nl-1),der,'.r')% red is new
hold on
plot(tp(1:nlp-1),derp,'.k')% black is old
xlabel('time, day')
ylabel('Time derivative of relative luminosity, 1/day')
title('Simulated transit lightcurve derivative')
print -deps derlumin.eps




%pause


% to make difference, need to find the index of t0 ~ 0



%figure(3)
%hold off
%plot(t,lum-lump,'.k')

%figure(4)
%hold off
%plot(t,der-derp,'.k')

%pause

save occsim2.mat