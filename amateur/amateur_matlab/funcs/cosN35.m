function ca = cosN35( v, w )
% Gives cosine of angle between each component of two vectors v(N,3,5) and w(N,3,5)
% It is assumed that the vectors are nonzero.
% Take care of NaN values!

ca = zeros(size(v));
for k=1:5
        ca(:,k) = cosN3(v(:,:,k),w(:,:,k));
        %ca(:,k) = (v(:,1,k).*w(:,1,k)+v(:,2,k).*w(:,2,k)+v(:,3,k).*w(:,3,k))./(norms(v).*norms(w));
end