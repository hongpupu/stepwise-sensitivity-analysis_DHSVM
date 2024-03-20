for a=1:7306 %loop time step (20 years-1st year for warming up)

%two searching curves
data = load(['C:\Users\WangLu\Desktop\Code\2 EFAST\dailyflow_example\daily.txt',num2str(a),'.txt']);
Data(1:6168,1) = data(1:6168,1);      
Data(1:6168,2) = data(6169:12336,1); 

N = 257;
NR = 2;  %Calculate each curve separately to obtain AV, AVi, and AVc, and then average the values of the two curves
  k=28;
  count = 1;
  OMi = 32;
  MI = 4;

  AV = 0;
  AVi = 0;
  AVci = 0;
  SI=zeros();
  TOTAL_SI=zeros();

for i=1:k %loop through parameters
    % Initialize AV,AVi,AVci to zero. 
    AV = 0;
    AVi = 0;
    AVci = 0;
    for L=1:NR%length(Y(1,t,u,i,:)) 
        Y1 = Data(count:count+N-1,L); 
        Y=Y1-mean(Y1);
        % Fourier coeff. at [1:OMi/2].
        NQ = (N-1)/2;
        N0 = NQ+1;
        COMPL = 0;
        Y_VECP = Y(N0+(1:NQ))+Y(N0-(1:NQ));
        Y_VECM = Y(N0+(1:NQ))-Y(N0-(1:NQ));
        for j=1:OMi/2
            ANGLE = j*2*(1:NQ)*pi/N;
            C_VEC = cos(ANGLE);
            S_VEC = sin(ANGLE);
            AC(j) = (Y(N0)+Y_VECP'*C_VEC')/N;
            BC(j) = Y_VECM'*S_VEC'/N;
            COMPL = COMPL+AC(j)^2+BC(j)^2;
        end
        % Computation of V_{(ci)}.
        Vci = 2*COMPL;
        AVci = AVci+Vci;
        % Fourier coeff. at [P*OMi, for P=1:MI].
        COMPL = 0;
        Y_VECP = Y(N0+(1:NQ))+Y(N0-(1:NQ));
        Y_VECM = Y(N0+(1:NQ))-Y(N0-(1:NQ));
        for j=OMi:OMi:OMi*MI
            ANGLE = j*2*(1:NQ)*pi/N;
            C_VEC = cos(ANGLE');
            S_VEC = sin(ANGLE');
            AC(j) = (Y(N0)+Y_VECP'*C_VEC)/N;
            BC(j) = Y_VECM'*S_VEC/N;
            COMPL = COMPL+AC(j)^2+BC(j)^2;
        end
        % Computation of V_i.
        Vi = 2*COMPL; %ith first-order variance
        AVi = AVi+Vi; %total first-order variance
        % Computation of the total variance
        % in the time domain.
        AV = AV+Y'*Y/N; %the total variance of the two search curves
    end %L
    % Computation of sensitivity indexes.
    %AV = AV/length(Y(1,t,1,i,:)); %AV=average total variance
    %AVi = AVi/length(Y(1,t,1,i,:));
    %AVci = AVci/length(Y(1,t,1,i,:));
    
      count = count + N;
    AV = AV/NR;        %the average of the two curves
    AVi = AVi/NR;
    AVci = AVci/NR;
    Si(i) = AVi/AV;
    Sti(i) = 1-AVci/AV;
end %i

   xlswrite('C:\Users\WangLu\Desktop\sensitivity_index.xlsx',Si);

end %each time step