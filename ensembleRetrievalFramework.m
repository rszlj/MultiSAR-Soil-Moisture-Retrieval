function [mvR,sig,X, cov_ave] = ensembleRetrievalFramework(bso,ensembleMode,forwardCube_IIEM,Ne,Nc,Nt,feMode,dryDown)
%UNTITLED Summary of this function goes here
%   
%   bso: input backscatter time series
%   ensembleMode:
%           1: snapshot ensemble
%           2: muli-temporal ensemble with full series
%           3: muli-temporal ensemble with random temporal sampling
%   forwardCube_IIEM: forward LUTs
%   Ne: the number of sub-retrievals
%   Nc: the number of selected channels  
%   Nt: the number of selected time instances
%   feMode: the mode of feature selction, 
%           1: randomly select Nc channels and/or Nt time instances, 
%           0: randonly select 1-Nc channels and/or 1-Nt time instances
%   dryDown:
%           1: dry down constrant is applied
%           0: dry down constrant is not applied


rng(88,'twister');

dateNum=length(unique(bso(:,4)));
mvR=cell(1,Ne);
sig=nan(1,Ne);
X=nan(1,Ne);

if ensembleMode==1 % snapshot ensemble retreival
    [tindex] = randomFeatures(bso,Nc,Ne,feMode);% randomly select channels
    for i=1:Ne
        tsig=nan(1,dateNum); %initalize the rms height series
        tmv=tsig; %initalize the mv series
        tX=tsig; %initalize the correlation length or VWC series 
        for k=1:dateNum % retrieve soil moisture day by day
            index4dayk=((bso(:,4)==k)+tindex{i})==2;
            bso4dayk=bso(index4dayk,:); % the selected channles for day k
            [ tsig(k),tmv(k),tX(k),~,~ ] = mvRetrieval_3in1('bare',Nt,bso4dayk,forwardCube_IIEM, dryDown); %the number of date is 1 for snapshot
        end
        mvR{i}=tmv;
        sig(i)=mean(tsig);
        X(i)=mean(tX);
    end
elseif ensembleMode==2 % multi-temporal with full series 
    dateNum=length(unique(bso(:,4))); % length of time series
    [tindex] = randomFeatures(bso,Nc,Ne,feMode);
    for i=1:Ne
        bso4i=bso(tindex{i}==1,:);
        [ tsig,mvR{i},tX,~,~ ] = mvRetrieval_3in1('bare',dateNum,bso4i,forwardCube_IIEM, dryDown);
        sig(i)=tsig(1); % tsig is a time series rms height with a same value because of the time-invariant roughness assumption
        X(i)=tX(1); % tX is a time series correlation length or vegetation water conent with a same value because of the time-invariant roughness or vegetation assumption
    end
elseif ensembleMode==3
    dateNum=length(unique(bso(:,4)));
    [temporalIndex,DOYLUT]=temporalSampling(bso,Nt,Ne,feMode);
    [feIndex] = randomFeatures(bso,Nc,Ne,feMode);
    for i=1:Ne
        bso4i=bso((temporalIndex{i}+feIndex{i})==2,:);
        [ bso4i ] = reEncode( DOYLUT{i},bso4i );
        selectedDateNum=length(unique(bso4i(:,4))); % the number of real selected date,selectedDateNum=Nt if feMode = 1
        mv4fullSeries=nan(1,dateNum);
        [ tsig,tmv,tX,~,~ ] = mvRetrieval_3in1('bare',selectedDateNum,bso4i,forwardCube_IIEM, dryDown);
        mv4fullSeries(DOYLUT{i}(1,:))=tmv;% since Nt~=dateNum, soil moisture retrieval is only made on the selected dates, with the unselected being NaN. 
        mvR{i}=mv4fullSeries;
        sig(i)=tsig(1);
        X(i)=tX(1);
    end
else
    fprintf('please select the ensemble mode: \n1: snapshot merge; \n2: multi-temporal with full series; \n3: muli-temporal ensemble with random temporal sampling\n')
end
[cov_ave]=pairWise_covariance(mvR);% calculate averaged pair wise covariance
[mvR] = mergeMethod(mvR,1);% merge the sub-retrievals
end

function [index,DOYLUT]=temporalSampling(bso,Nt,Ne,feMode)
%generate the random dates for ensemble
dateN=unique(bso(:,4));
dateNum=length(dateN);
DOYLUT=cell(Ne,1);
index=cell(Ne,1);

if Nt~=dateNum
   
tindex=zeros(Ne,dateNum);
    flag=1;
    while flag==1
        if feMode==1
            tempN1=1+round(rand(1,Ne)*(Nt-1));
        else
            tempN1=ones(1,Ne)*Nt;
        end
        tempN2=rand(Ne,dateNum);
        [~,id]=sort(tempN2,2,'descend');   
        for i=1:Ne
        tindex(i,id(i,1:tempN1(i)))=1;
        DOYLUT{i}=[dateN(sort(id(i,1:tempN1(i))))';1:tempN1(i)];
        tempppp=zeros(length(bso(:,4)),1);
            for j=1:length(DOYLUT{i}(1,:))
            tempppp(bso(:,4)==DOYLUT{i}(1,j))=1;
            end
        index{i}=tempppp;    
        end
        flag=(sum(sum(tindex)==0)~=0);
    end
else
    for i=1:Ne
    DOYLUT{i}=[1:dateNum;1:dateNum];
    index{i}=ones(size(bso(:,4)));
    end
end

end
function [selectedFe] = randomFeatures(bso,Nc,Ne,femode)
%randomly select channles for each date of Ne sub-retrievals
%   selectedFe: 
dateN=unique(bso(:,4));
selectedFe=cell(Ne,1);
for i=1:Ne
    selectedFe4i=zeros(length(bso(:,1)),1);
    for k=1:length(dateN)
        tempid=find(bso(:,4)==dateN(k));
        polN=length(tempid);% the available channels of date k
        tpolMaxN=min(polN,Nc);% set Nc to the polN if Nc is larger than the number of available channels
        if femode==1
            tempN1=1+round(rand(1,1)*(tpolMaxN-1));
            tempN2=rand(1,polN);
        else
            tempN1=tpolMaxN;
            tempN2=rand(1,polN);
        end
        [~,id]=sort(tempN2,'descend');
        tempid=tempid(id(1:tempN1),1);
        selectedFe4i(tempid)=1;
    end
    selectedFe{i}=selectedFe4i;
end
end
function [ bso ] = reEncode( doyLUT, bso )
%Link between the doy system and temporal id
%   Detailed explanation goes here
doy4transfer=bso(:,4);
id1=unique(doy4transfer);
tid=nan(size(doy4transfer));
for i=1:length(id1)
   index1=doyLUT(1,:)==id1(i);
   index2=doy4transfer==id1(i);
   tid(index2)=doyLUT(2,index1);
end
bso(:,4)=tid;
end
function [mvRMerge]=mergeMethod(mvR,mode)
% merge the sub-retrievals

mvR=mvR(:);
mvR=cell2mat(mvR);

if mode==1%average
    mvRMerge=nanmean(mvR);   
elseif mode==2%medium
    mvRMerge=nanmedian(mvR); 
elseif mode==3
    tempstd=std(mvR);
    tempmv=mean(mvR);
    index1=mvR>(tempmv-tempstd);
    index2=mvR<(tempmv+tempstd);
    mvR((index1+index2)~=2)=nan;
    mvRMerge=nanmean(mvR);
end

end
function [cov_ave]=pairWise_covariance(mvR)
[simuNum,Ne]=size(mvR);
temp_corr = nan(Ne*(Ne-1),1);
k = 1;
for i = 1: Ne
    for j = 1:Ne
        if i ~=j
            temp_mv1 = cell2mat(mvR(:,i));
            temp_mv2 = cell2mat(mvR(:,j));
            temp = cov(temp_mv1(:),temp_mv2(:));
            temp_corr(k)= temp(1,2);
            k=k+1;
        end
    end
end
cov_ave = mean(temp_corr);
end
