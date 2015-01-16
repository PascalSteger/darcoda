%get number of total (sub/host)halos
%o=[];
%for i=1:length(hostno)
%    o(i)=length(hostno(subno==i-1));
%end

% exclude auto-correlation entries
bnzheqs=hostno==subno;
bnzhles=hostno<=subno;
bnzhges=hostno>=subno;
ops=[]; pos=[]; count=1;
for i=1:length(bnzheqs)-1
    pos(i)=0;
    ops(count)=0;
    if i==1 | hostno(i+1)>hostno(i)
         count=count+1;
    end
end
%get count of same entries in hostno
count=1;
for i=1:length(bnzheqs)-1
    ops(count)=ops(count)+1;
    if i==1 | hostno(i+1)>hostno(i)
        pos(i)=1;
        if i>1
            count=count+1;
        end
    end
end
opsacc=[]; ra=0;
for i=1:length(ops)
    opsacc(i)=ra;
    ra=ra+ops(i);
end

%find only the father (=next upper) host subhalo
