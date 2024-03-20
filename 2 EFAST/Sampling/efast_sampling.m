clear;
close all;
%reference:A quantitative model-independent method for global sensitivity
%analysis of model output
%% INPUT
NR = 2; %: no. of search curves_RESAMPLING 
k = 24; % # of input factors (parameters varied) 
NS = 257; % # of samples per search curve
wantedN=NS*k*NR; % wanted no. of sample points// total number of model evaluation

MI = 4; %: maximum number of fourier coefficients//文中的interference factor(usually 4 or higher)
% that may be retained in calculating the partial
% variances without interferences between the
% assigned frequencies

%% PARAMETERS AND ODE SETTINGS
Parameter_settings_EFAST;
% Computation of the frequency for the group
% of interest OMi and the # of sample points NS (here N=NS)
OMi = floor(((wantedN/NR)-1)/(2*MI)/k); %the largest among the set of wi frequencies;y = floor(x) 函数将x中元素取整,值y为不大于本身的最大整数。
NS = 2*MI*OMi+1;
if(NS*NR < 65)
    fprintf(['Error: sample size must be >= ' ...
    '65 per factor.\n']);
    return;
end
%% PARAMATERS DISTRIBUTION
for i=1:k
     OMci = SETFREQ(k,OMi/2/MI,i);
     for L=1:NR
        % Setting the vector of frequencies OM
        % for the k parameters
        cj = 1;
        for j=1:k
            if(j==i)
                % For the parameter (factor) of interest
                OM(i) = OMi;
            else
                % For the complementary group.
                OM(j) = OMci(cj);
                cj = cj+1;
            end
        end
        % Setting the relation between the scalar
        % variable S and the coordinates
        % {X(1),X(2),...X(k)} of each sample point.
        FI = rand(1,k)*2*pi; % random phase shift
        S_VEC = pi*(2*(1:NS)-NS-1)/NS;
        OM_VEC = OM(1:k);
        FI_MAT = FI(ones(NS,1),1:k)';
        ANGLE = OM_VEC'*S_VEC+FI_MAT;
        
        X(:,:,i,L) = 0.5+asin(sin(ANGLE'))/pi; % between 0 and 1
        
        % Transform distributions from standard
        % uniform to general.
        X(:,:,i,L) = parameterdist(X(:,:,i,L),pmax,pmin,0,1,NS,'unif');
        
     end
end
Y=squeeze(X);
