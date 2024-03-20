clear all
close all
path(path,'..\..\..\FUZZCLUST')
%data set
load percentage.txt
data.X = percentage(:,:);


%parameters
param.c=4;
param.m=2;
param.e=1e-6;
param.ro=ones(1,param.c);
param.val=3;
%normalization
data=clust_normalize(data,'range');
%clustering
result = GKclust(data,param);
plot(data.X(:,1),data.X(:,2),'b.',result.cluster.v(:,1),result.cluster.v(:,2),'ro');
hold on
%draw contour-map
new.X=data.X;
eval=clusteval(new,result,param);
%validation
result = validity(result,data,param);
result.validity
%PCA projection
param.q = 2;
result = PCA(data,param,result);
figure(2)
plot(result.proj.P(:,1),result.proj.P(:,2),'.')
hold on
plot(result.proj.vp(:,1),result.proj.vp(:,2),'r*')
perf = projeval(result,param);

perf