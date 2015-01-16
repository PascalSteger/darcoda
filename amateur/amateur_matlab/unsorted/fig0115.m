tidalfieldinput()

X=pos(1,:); Y=pos(2,:); Z=pos(3,:);
U1=ev(1,:); V1=ev(2,:); W1=ev(3,:);

T1=tens(1,:); T2=tens(2,:); T3=tens(3,:);
T4=tens(4,:); T5=tens(5,:); T6=tens(6,:);

scatter3(X,Y,Z,ovdens/1000,'.')
scatter3(X,Y,Z,dens/100,'.')
hold on
quiver3(X,Y,Z,U1,V1,W1)
hold off