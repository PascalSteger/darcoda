function ns = norma( A )
% ns(i) gives length of subvectors A(i,:)
% where A has dimensions (:,3)

%ns = zeros(length(A(:,1)),1);
%for i=1:length(A(:,1))
%        ns(i) = norm(A(i,:));
%end

ns = normv( A(:,1), A(:,2), A(:,3) );
