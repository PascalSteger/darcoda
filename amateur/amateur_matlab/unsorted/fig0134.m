scatter(hp0cmx,hp0cmy,'r.');
hold on; 
%quiver(hp0cmx,hp0cmy,hp0vcmx,hp0vcmy,'b');
%quiver(hp0cmx,hp0cmy,tt(1,:,1),tt(2,:,1),'r');
quiver(hp0cmx,hp0cmy,hp_l(1,:,1)./hp_ltot(1,:,1),hp_l(2,:,1)./hp_ltot(1,:,1),'m');
quiver(hp0cmx,hp0cmy,hp_l(1,:,2)./hp_ltot(1,:,2),hp_l(2,:,2)./hp_ltot(1,:,2),'r');
quiver(hp0cmx,hp0cmy,hp0vcmx,hp0vcmy,'g');

axis square; hold off;