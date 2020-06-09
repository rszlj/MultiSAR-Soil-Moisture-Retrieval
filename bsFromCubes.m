function [ thh,thv,tvv] = bsFromCubes(sig,cl,mv,thetai,cubeL )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
interpN=length(sig);
thh=nan(interpN,1);
thv=nan(interpN,1);
tvv=nan(interpN,1);

paraS=cl;
paraSN=cubeL.XN;
sigN=cubeL.sigN;
mvN=cubeL.mvN;


incid=cubeL.incid;
for i=1:interpN
    temp1=find(sigN>=sig(i));
    temp2=find(paraSN>=paraS(i));
    temp3=find(mvN>=mv(i));
    id=temp1(1)+(temp2(1)-1)*64+(temp3(1)-1)*64*64;
    temphh=cubeL.hh(:,id);
    temphv=cubeL.hv(:,id);
    tempvv=cubeL.vv(:,id);
    thh(i)=interp1(incid',temphh,thetai,'spline');
    thv(i)=interp1(incid',temphv,thetai,'spline');
    tvv(i)=interp1(incid',tempvv,thetai,'spline');
end
end
