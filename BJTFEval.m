%% Box-Jenkins model evaluation only for TF part
%look up table generation
function  M=BJTFEval(model,u,y)
    
    if strcmp(class(model),'idtf')==1
        numBJ=model.Numerator;
        denBJ=model.Denominator;
    else
    %box jenkins model M3 just with tf part consideration
        numBJ=model.B;
        denBJ=model.F;
    end
   
    %evaluating the difference equation
    yMean=mean(y); uMean=mean(u);
    %fourth order BJ transfer function
    uk=u-uMean;                           %input vector minus mean
    yk=y-yMean;                           %output vector minus mean
    ykBJ=zeros(1,length(yk));       
    ykBJ(1:5)=yk(1:5);                %BJ just tf output vector first 5 elements equal
    ykBJTable=ykBJ;                 %Look up table with unmodeled behavior of the model

    for i=6:length(y)

        ykA= numBJ(1)*uk(i)+ numBJ(2)*uk(i-1) + numBJ(3)*uk(i-2) +numBJ(4)*uk(i-3)+numBJ(5)*uk(i-4); %input based BJ model
        ykB= denBJ(2)*yk(i-1) + denBJ(3)*yk(i-2) +denBJ(4)*yk(i-3)+denBJ(5)*yk(i-4);                 %output based BJ model
        ykBJ(i)=ykA-ykB;

        ykBJTable(i)=y(i)-ykBJ(i)-yMean;  %look up table calculation

    end
    ykBJTable(1:5)=y(1:5);                %BJ just tf output vector first 5 elements equal
    M=round(ykBJTable',1);

end