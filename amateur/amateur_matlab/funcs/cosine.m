function ca = cosine( vx, vy, vz, wx, wy, wz )
% Gives cosine of angle between two vectors (vx,vy,vz) and (wx,wy,wz)
% It is assumed that the vectors are nonzero.
% Take care of NaN values!

ca = (vx.*wx + vy.*wy + vz.*wz)./(normN1(vx,vy,vz).*normN1(wx,wy,wz));
