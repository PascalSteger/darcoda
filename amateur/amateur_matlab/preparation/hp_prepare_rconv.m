%amautinput()
%ahfinput() %for sub->host matching

cmx=[]; cmy=[]; cmz=[]; cmxold=[]; cmyold=[]; cmzold=[];
for k=1:length(hp0cm)
    cmx(k)=(hp0cm(1,k)-hp0cm(1,hostno(k)+1));
     cmy(k)=(hp0cm(2,k)-hp0cm(2,hostno(k)+1));
      cmz(k)=(hp0cm(3,k)-hp0cm(3,hostno(k)+1));
    cmxold(k)=(hp0cm(1,k)-hp0cm(1,1));
     cmyold(k)=(hp0cm(2,k)-hp0cm(2,1));
      cmzold(k)=(hp0cm(3,k)-hp0cm(3,1));
end

hp0jx=hp0j(1,:);
hp0jy=hp0j(2,:);
hp0jz=hp0j(3,:);
hp1jx=hp1j(1,:);
hp1jy=hp1j(2,:);
hp1jz=hp1j(3,:);
hp2jx=hp2j(1,:);
hp2jy=hp2j(2,:);
hp2jz=hp2j(3,:);
hp3jx=hp3j(1,:);
hp3jy=hp3j(2,:);
hp3jz=hp3j(3,:);
hp4jx=hp4j(1,:);
hp4jy=hp4j(2,:);
hp4jz=hp4j(3,:);
hp5jx=hp5j(1,:);
hp5jy=hp5j(2,:);
hp5jz=hp5j(3,:);

hp0vx=hp0evmax(1,:);
hp0vy=hp0evmax(2,:);
hp0vz=hp0evmax(3,:);
hp1vx=hp1evmax(1,:);
hp1vy=hp1evmax(2,:);
hp1vz=hp1evmax(3,:);
hp2vx=hp2evmax(1,:);
hp2vy=hp2evmax(2,:);
hp2vz=hp2evmax(3,:);
hp3vx=hp3evmax(1,:);
hp3vy=hp3evmax(2,:);
hp3vz=hp3evmax(3,:);
hp4vx=hp4evmax(1,:);
hp4vy=hp4evmax(2,:);
hp4vz=hp4evmax(3,:);
hp5vx=hp5evmax(1,:);
hp5vy=hp5evmax(2,:);
hp5vz=hp5evmax(3,:);

hp0jtot=sqrt(hp0jx.*hp0jx+hp0jy.*hp0jy+hp0jz.*hp0jz);
hp1jtot=sqrt(hp1jx.*hp1jx+hp1jy.*hp1jy+hp1jz.*hp1jz);
hp2jtot=sqrt(hp2jx.*hp2jx+hp2jy.*hp2jy+hp2jz.*hp2jz);
hp3jtot=sqrt(hp3jx.*hp3jx+hp3jy.*hp3jy+hp3jz.*hp3jz);
hp4jtot=sqrt(hp4jx.*hp4jx+hp4jy.*hp4jy+hp4jz.*hp4jz);
hp5jtot=sqrt(hp5jx.*hp5jx+hp5jy.*hp5jy+hp5jz.*hp5jz);

hp0vtot=sqrt(hp0vx.*hp0vx+hp0vy.*hp0vy+hp0vz.*hp0vz);
hp1vtot=sqrt(hp1vx.*hp1vx+hp1vy.*hp1vy+hp1vz.*hp1vz);
hp2vtot=sqrt(hp2vx.*hp2vx+hp2vy.*hp2vy+hp2vz.*hp2vz);
hp3vtot=sqrt(hp3vx.*hp3vx+hp3vy.*hp3vy+hp3vz.*hp3vz);
hp4vtot=sqrt(hp4vx.*hp4vx+hp4vy.*hp4vy+hp4vz.*hp4vz);
hp5vtot=sqrt(hp5vx.*hp5vx+hp5vy.*hp5vy+hp5vz.*hp5vz);

rtot=sqrt(cmx.*cmx+cmy.*cmy+cmz.*cmz);
rtotold=sqrt(cmxold.*cmxold+cmyold.*cmyold+cmzold.*cmzold);

%autocorrelations of angular momenta
cos01=(hp0jx.*hp1jx+hp0jy.*hp1jy+hp0jz.*hp1jz)./(hp0jtot.*hp1jtot);
cos02=(hp0jx.*hp2jx+hp0jy.*hp2jy+hp0jz.*hp2jz)./(hp0jtot.*hp2jtot);
cos03=(hp0jx.*hp3jx+hp0jy.*hp3jy+hp0jz.*hp3jz)./(hp0jtot.*hp3jtot);
cos04=(hp0jx.*hp4jx+hp0jy.*hp4jy+hp0jz.*hp4jz)./(hp0jtot.*hp4jtot);
cos05=(hp0jx.*hp5jx+hp0jy.*hp5jy+hp0jz.*hp5jz)./(hp0jtot.*hp5jtot);

cos12=(hp1jx.*hp2jx+hp1jy.*hp2jy+hp1jz.*hp2jz)./(hp1jtot.*hp2jtot);
cos13=(hp1jx.*hp3jx+hp1jy.*hp3jy+hp1jz.*hp3jz)./(hp1jtot.*hp3jtot);
cos14=(hp1jx.*hp4jx+hp1jy.*hp4jy+hp1jz.*hp4jz)./(hp1jtot.*hp4jtot);
cos15=(hp1jx.*hp5jx+hp1jy.*hp5jy+hp1jz.*hp5jz)./(hp1jtot.*hp5jtot);

cos23=(hp2jx.*hp3jx+hp2jy.*hp3jy+hp2jz.*hp3jz)./(hp2jtot.*hp3jtot);
cos24=(hp2jx.*hp4jx+hp2jy.*hp4jy+hp2jz.*hp4jz)./(hp2jtot.*hp4jtot);
cos25=(hp2jx.*hp5jx+hp2jy.*hp5jy+hp2jz.*hp5jz)./(hp2jtot.*hp5jtot);

cos34=(hp3jx.*hp4jx+hp3jy.*hp4jy+hp3jz.*hp4jz)./(hp3jtot.*hp4jtot);
cos35=(hp3jx.*hp5jx+hp3jy.*hp5jy+hp3jz.*hp5jz)./(hp3jtot.*hp5jtot);

cos45=(hp4jx.*hp5jx+hp4jy.*hp5jy+hp4jz.*hp5jz)./(hp4jtot.*hp5jtot);

% cos(angle between position vector from subhalos to host halo and angular momentum vector)
cosrj0=(cmx.*hp0jx+cmy.*hp0jy+cmz.*hp0jz)./(hp0jtot.*rtot);
cosrj1=(cmx.*hp1jx+cmy.*hp1jy+cmz.*hp1jz)./(hp1jtot.*rtot);
cosrj2=(cmx.*hp2jx+cmy.*hp2jy+cmz.*hp2jz)./(hp2jtot.*rtot);
cosrj3=(cmx.*hp3jx+cmy.*hp3jy+cmz.*hp3jz)./(hp3jtot.*rtot);
cosrj4=(cmx.*hp4jx+cmy.*hp4jy+cmz.*hp4jz)./(hp4jtot.*rtot);
cosrj5=(cmx.*hp5jx+cmy.*hp5jy+cmz.*hp5jz)./(hp5jtot.*rtot);

%the same, but from subhalo to big host halo only
cosrj0old=(cmxold.*hp0jx+cmyold.*hp0jy+cmzold.*hp0jz)./(hp0jtot.*rtotold);
cosrj1old=(cmxold.*hp1jx+cmyold.*hp1jy+cmzold.*hp1jz)./(hp1jtot.*rtotold);
cosrj2old=(cmxold.*hp2jx+cmyold.*hp2jy+cmzold.*hp2jz)./(hp2jtot.*rtotold);
cosrj3old=(cmxold.*hp3jx+cmyold.*hp3jy+cmzold.*hp3jz)./(hp3jtot.*rtotold);
cosrj4old=(cmxold.*hp4jx+cmyold.*hp4jy+cmzold.*hp4jz)./(hp4jtot.*rtotold);
cosrj5old=(cmxold.*hp5jx+cmyold.*hp5jy+cmzold.*hp5jz)./(hp5jtot.*rtotold);

% cos(angle between position vector from biggest halo to smaller halos and major axis)
cosrv0=(cmx.*hp0vx+cmy.*hp0vy+cmz.*hp0vz)./(hp0vtot.*rtot);
cosrv1=(cmx.*hp1vx+cmy.*hp1vy+cmz.*hp1vz)./(hp1vtot.*rtot);
cosrv2=(cmx.*hp2vx+cmy.*hp2vy+cmz.*hp2vz)./(hp2vtot.*rtot);
cosrv3=(cmx.*hp3vx+cmy.*hp3vy+cmz.*hp3vz)./(hp3vtot.*rtot);
cosrv4=(cmx.*hp4vx+cmy.*hp4vy+cmz.*hp4vz)./(hp4vtot.*rtot);
cosrv5=(cmx.*hp5vx+cmy.*hp5vy+cmz.*hp5vz)./(hp5vtot.*rtot);

% cos(angle between position vector from biggest halo to smaller halos and major axis)
cosrv0old=(cmxold.*hp0vx+cmyold.*hp0vy+cmzold.*hp0vz)./(hp0vtot.*rtotold);
cosrv1old=(cmxold.*hp1vx+cmyold.*hp1vy+cmzold.*hp1vz)./(hp1vtot.*rtotold);
cosrv2old=(cmxold.*hp2vx+cmyold.*hp2vy+cmzold.*hp2vz)./(hp2vtot.*rtotold);
cosrv3old=(cmxold.*hp3vx+cmyold.*hp3vy+cmzold.*hp3vz)./(hp3vtot.*rtotold);
cosrv4old=(cmxold.*hp4vx+cmyold.*hp4vy+cmzold.*hp4vz)./(hp4vtot.*rtotold);
cosrv5old=(cmxold.*hp5vx+cmyold.*hp5vy+cmzold.*hp5vz)./(hp5vtot.*rtotold);
