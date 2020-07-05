%% Box-Jenkins model evaluation only for TF part
function  M=BJTFEvalRand(model,u,y,nRounded)
    B=model.B;
    F=model.F;
    BD=B;
    FD=F;
   %evaluating the difference equation
    yMean=mean(y); uMean=mean(u);
    %fourth order BJ transfer function
    uk=u;                           %input vector minus mean
    yk=y;                     
    ykBJ=zeros(1,length(yk));       
    ykBJ(1:length(FD))=yk(1:length(FD));                %BJ just tf output vector first 5 elements equal
    ykBJTable=ykBJ;                 %Look up table with unmodeled behavior of the model

    for i=length(FD)+1:length(y)
        ykU=0; ykE=0;ykY=0;
        
        for j=0:length(BD)-1
            ykU=BD(j+1)*uk(i-j)  +ykU; %input
             ykE=CF(j+1)*e(i-j)   +ykE; %white noise
        end
        for k=2:length(FD)
             ykYa=-FD(k)*yk(i-k+1);    %output
             ykY=ykYa+ykY;
             [i k i-k+1 ykY];
        end
        ykBJ(i)=ykU+ykE+ykY;
  
    end
    ykBJTable=y-ykBJ';%-yMean;  %look up table calculation
    ykBJTable(1:length(FD),1)=y(1:length(FD),1);                %BJ just tf output vector first 5 elements equal
    ykBJTable=ykBJTable';
    M=round(ykBJTable',nRounded);
end
