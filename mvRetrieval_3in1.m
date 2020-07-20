function [ sig,mv,X,resError,resEmat ] = mvRetrieval_3in1(lcType,dateNum,bso,forwardCube_IIEM, dryDown)
%Soil moisture retrieval methods 3 in 1
%   lcType: landcover type, i.e. 1 bare; 2 wheat; 3 grass
%   bso: backscattering coefficient series
%   imageNum: the number of dates, each date may have multiple acquisitions
%   lia: incidence angle for each backscattering coefficient
%   bsID: the id of polarization, 1 HH; 2 VH; 3 VV;
%   tID: date id for each polarzation
%   frID: the frequency id for each polarzation, 1 L; 2 C; 3 X 



%Parameter initialization
NVAR_sig = 1;           % Number of unknowns describing the rms height, being 1 for time-invariant roughness
NVAR_X = 1;           % Number of unknowns describing the VWC or correlation length, being 1 for time-invariant roughness or vegetation
NVAR_mv = dateNum;           % Number of soil moisture unknowns, being same to the number of dates
NVAR=NVAR_sig+NVAR_X+NVAR_mv;

NIND = 20;           % Number of individuals per subpopulations. Not that this is not Ne 
MAXGEN = 50;        % maximum Number of generations
GGAP = .9;           % Generation gap, how many new individuals are created
PRECI = 8;          % Precision of binary representation

% initialization unknowns,you may indicate the range of each unknonw here
if strcmp(lcType,'wheat')
    FieldD_X = [rep([PRECI],[1, NVAR_X]); rep([0.5;4],[1, NVAR_X]);...
          rep([1; 0; 1 ;1], [1, NVAR_X])]; % range of VWC is 0.5 to 4 kg/m2
    FieldD_sig = [rep([PRECI],[1, NVAR_sig]); rep([0.5;4],[1, NVAR_sig]);...
          rep([1; 0; 1 ;1], [1, NVAR_sig])];% range of sig is 0.5 to 4 cm
    FieldD_mv = [rep([PRECI],[1, NVAR_mv]); rep([0.04;0.42],[1, NVAR_mv]);...
          rep([1; 0; 1 ;1], [1, NVAR_mv])]; % range of mv is 0.04 to 0.42 m3/m3
elseif strcmp(lcType,'grass')
    FieldD_X = [rep([PRECI],[1, NVAR_X]); rep([0.1;2],[1, NVAR_X]);...
          rep([1; 0; 1 ;1], [1, NVAR_X])];  % range of VWC is 0.1 to 2 kg/m2
    FieldD_sig = [rep([PRECI],[1, NVAR_sig]); rep([0.5;3],[1, NVAR_sig]);...
          rep([1; 0; 1 ;1], [1, NVAR_sig])];% range of sig is 0.5 to 3 cm
    FieldD_mv = [rep([PRECI],[1, NVAR_mv]); rep([0.04;0.42],[1, NVAR_mv]);...
          rep([1; 0; 1 ;1], [1, NVAR_mv])]; % range of mv is 0.04 to 0.42 m3/m3
else
    FieldD_X = [rep([PRECI],[1, NVAR_X]); rep([5;35],[1, NVAR_X]);...
          rep([1; 0; 1 ;1], [1, NVAR_X])]; % range of correlation length is 5 to 35 cm
    FieldD_sig = [rep([PRECI],[1, NVAR_sig]); rep([0.5;3],[1, NVAR_sig]);...
          rep([1; 0; 1 ;1], [1, NVAR_sig])];% range of correlation length is 0.5 to 3 cm
    FieldD_mv = [rep([PRECI],[1, NVAR_mv]); rep([0.04;0.45],[1, NVAR_mv]);...
          rep([1; 0; 1 ;1], [1, NVAR_mv])]; % range of mv is 0.04 to 0.45 m3/m3
end

% Initialise population
   Chrom = crtbp(NIND, NVAR*PRECI);
   
% Reset counters
   Best = NaN*ones(MAXGEN,1);	% best in current population
   gen = 0;			% generational counter

% Decoding the Chrom to real value matrix
 [ sig,X,mv ] = bs2rv_2n1( Chrom,FieldD_sig,FieldD_X,FieldD_mv);
 
% Evaluate initial population
   [ ObjV,resEmat ] = Obj4multiFrequency( sig,X,mv,bso,forwardCube_IIEM, dryDown );
        
objectV=nan(MAXGEN,1);
k=1;
% Generational loop
   while gen < MAXGEN

    % Assign fitness-value to entire population
       FitnV = ranking(ObjV);

    % Select individuals for breeding
       SelCh = select('sus', Chrom, FitnV, GGAP);

    % Recombine selected individuals (crossover)
       SelCh = recombin('xovsp',SelCh,0.7);

    % Perform mutation on offspring
       SelCh = mut(SelCh);

    % Decoding the Chrom to real value matrix and Evaluate offspring, call objective function
    [ sig,X,mv ] = bs2rv_2n1( SelCh,FieldD_sig,FieldD_X,FieldD_mv);
    [ ObjVSel,resEmat ] = Obj4multiFrequency( sig,X,mv,bso,forwardCube_IIEM, dryDown );

    objectV(k)=min(ObjVSel);
    k=k+1;
    % Reinsert offspring into current population
       [Chrom, ObjV]=reins(Chrom,SelCh,1,1,ObjV,ObjVSel);
    % Increment generational counter
       gen = gen+1;

   end 
% End of GA
   [ sig,X,mv ] = bs2rv_2n1( Chrom,FieldD_sig,FieldD_X,FieldD_mv);
   [ ObjV,resEmat ] = Obj4multiFrequency( sig,X,mv,bso,forwardCube_IIEM, dryDown);
   [id,value]=find(ObjV==min(ObjV));
   sig=sig(id(1),:);
   mv=mv(id(1),:);
   X=X(id(1),:);
   resError=min(ObjV);
   resEmat=resEmat(id(1),:);
end
function [ sig,cl,mv ] = bs2rv_2n1( Chrom,FieldD_sig,FieldD_cl,FieldD_mv)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%decoding sig
[~,sigNum]=size(FieldD_sig);
bN1=FieldD_sig(1)*sigNum;
sig=bs2rv(Chrom(:,1:bN1),FieldD_sig);

%decoding cl
[~,clNum]=size(FieldD_cl);
bN2=FieldD_cl(1)*clNum+bN1;
cl=bs2rv(Chrom(:,bN1+1:bN2),FieldD_cl);

%decoding mv
[~,mvNum]=size(FieldD_mv);
% bN2=Field_mv(1)*mvNum;
mv=bs2rv(Chrom(:,bN2+1:end),FieldD_mv);

varNum=[sigNum,clNum,mvNum];
maxValue=max(varNum);
if sigNum~=maxValue
    sig=repmat(sig,1,maxValue);
end
if clNum~=maxValue
    cl=repmat(cl,1,maxValue);
end
if mvNum~=maxValue
    mv=repmat(mv,1,maxValue);
end
end
function [ ObjV,resError ] = Obj4multiFrequency( sig,X,mv,bso,forwardCube_IIEM, dryDown )
%Cost function for multi-temporal retrieval
%   Detailed explanation goes here

%bs lia pol dateID sideID frequency
bs=bso(:,1);
lia=bso(:,2);
bsID=bso(:,3);
tID=bso(:,4);
sideID=bso(:,5);
frID=bso(:,6);

bsNum=length(bs);

[m,n]=size(mv);
unitID=unique(tID);
simBs=nan(m,bsNum);

allID=frID*1000+tID*10+sideID;
uniqueID=unique(allID);
uniqueIDnum=length(uniqueID);

uniquetID=unique(tID);
for i=1:uniqueIDnum
    index=find(allID==uniqueID(i));%+++++++++++++++++++++++++++++++++++++++++=
    tempfr=floor(uniqueID(i)/1000);
    temptID=mod(floor(uniqueID(i)/10),100);
    [temptID,~]=find(uniquetID==temptID);
    for j=1:m
        if tempfr==1
            cube=forwardCube_IIEM.cubeL;
        elseif tempfr==2
            cube=forwardCube_IIEM.cubeC;
        else
            cube=forwardCube_IIEM.cubeX;
        end
         [ tempbs ] = cubeInterpolation( sig(j,temptID),X(j,temptID),mv(j,temptID),lia(index(1)),cube,tempfr);
         simBs(j,index)=tempbs(bsID(index));
    end
end


bs=repmat(bs',m,1);
if dryDown==1
    [~,sq3]=sort(mv,2,'descend'); %Zhu et al., 2019
    mvFac=sum(abs(repmat(1:n,m,1)-sq3),2);
    mvFac(mvFac==0)=1;
else
    mvFac = 1;
end
ObjV=sqrt(nanmean((simBs-bs).^2,2)).*mvFac;
resError=abs(simBs-bs);
end
function [ bs ] = cubeInterpolation( sig,X,mv,thetai,cube,tempfr)
%A simple two 1-D interpolation functions to get the sigmanull from the LUT
%   Note that 2-D interpolation should be used in real applications
sigN=cube.sigN;
mvN=cube.mvN;
paraS=X;
paraSN=cube.XN;
incid=cube.incid;
interpN=length(sig);
thh=nan(interpN,1);
thv=nan(interpN,1);
tvv=nan(interpN,1);

for i=1:interpN
    temp1=find(sigN>=sig(i));
    temp2=find(paraSN>=paraS(i));
    temp3=find(mvN>=mv(i));
    id=temp1(1)+(temp2(1)-1)*64+(temp3(1)-1)*64*64;
    if tempfr==3
    temphh=cube.hh(:,id);
    thh(i)=interp1(incid',temphh,thetai,'spline'); 
    thv(i)=thh(i);
    tvv(i)=tvv(i);
    else
    temphh=cube.hh(:,id);
    temphv=cube.hv(:,id);
    tempvv=cube.vv(:,id);
    thh(i)=interp1(incid',temphh,thetai,'spline');
    thv(i)=interp1(incid',temphv,thetai,'spline');
    tvv(i)=interp1(incid',tempvv,thetai,'spline');
    end
end
bs=[thh,thv,tvv];
end