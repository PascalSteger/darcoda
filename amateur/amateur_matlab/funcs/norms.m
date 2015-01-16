function ns = norms( A )
% ns(:,1,j) gives length of subvectors A(:,1:3,j)
% where A has dimensions (N,3,5)

ns = zeros(length(A(:,1,1)),1,5);
for i=1:5
        for j=1:length(A(:,1,1));
                ns(j,1,i) = norm([A(j,1,i),A(j,2,i),A(j,3,i)]);
        end
end
