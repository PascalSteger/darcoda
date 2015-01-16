% fig0115.m: visualization of tidal field vectors, cf. Hahn 2007, fig. 8
clf
set(gca,'FontSize',15);

X=Tpos(1,:); Y=Tpos(2,:); Z=Tpos(3,:); %xyz,#halo
evax = Teva(1,:); evay = Teva(2,:); evaz = Teva(3,:);

T1=Ttens(1,:); T2=Ttens(2,:); T3=Ttens(3,:);
T4=Ttens(4,:); T5=Ttens(5,:); T6=Ttens(6,:);

% select some 2000 points to be displayed
n = size(Tpos(1,:));
n = n(2);
f = ceil(n.*rand(2000,1));


% and plot them, with size propto overdensity
subplot(1,2,1);
scatter3(X(f),Y(f),Z(f),Tovdens(f),'.')
hold on
% as well as the direction of the main eigenvector of the tidal field
% the direction should not enter here
quiver3(X(f),Y(f),Z(f),evax(f),evay(f),evaz(f),'-')
hold off

% now take a slice
subplot(1,2,2);
f = prctile(Z,[40 60]);
f = f(1)<Z & Z<f(2);
quiver(X(f),Y(f),evax(f),evay(f),'-'); hold on;
scatter(X(f),Y(f),Tovdens(f),'.');
hold off;