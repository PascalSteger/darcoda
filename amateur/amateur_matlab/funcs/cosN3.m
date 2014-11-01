function ca = cosN3( v, w )
% Gives cosine of angle between each component of two vectors v(N,3) and w(N,3)
% It is assumed that the vectors are nonzero.
% Take care of NaN values!

ca = cosine(v(:,1),v(:,2),v(:,3),w(:,1),w(:,2),w(:,3));
