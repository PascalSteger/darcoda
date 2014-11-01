function ca = vcos( v, w )
% Gives cosine of angle between two 3D-vectors v and w

ca = (v(1).*w(1) + v(2).*w(2) + v(3).*w(3))./(norm(v).*norm(w));