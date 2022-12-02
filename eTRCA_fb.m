function [acc] = eTRCA_fb(XTrain,YTrain,XTest,YTest,num_fb)
%% output the accuracy of eTRCA with filter bank
%   XTrain: num_train * num_channel * num_sampls * num_fb
%   YTrain: num_train * 1
%   XTest: num_test * num_channel * num_sampls * num_fb
%   YTest: num_test * 1
%   num_fb: number of the filter bank


[~, num_channel, num_sampls] = size(XTrain(:,:,:,1));
class_all = unique(YTrain);
num_class = length(class_all);
num_test=size(XTest,1);
crall=zeros(num_test,num_fb,num_class);
fb = (1:num_fb).^(-1.25)+0.25;
for fb_i=1:num_fb
    XTrain_fbi=XTrain(:,:,:,fb_i);
    XTest_fbi=XTest(:,:,:,fb_i);
    Template=cell(1,num_class);
    for class_i = 1:num_class
        XTrain_class_i = squeeze(XTrain_fbi(YTrain==class_all(class_i),:,:));
        Template{class_i} = squeeze(mean(XTrain_class_i,1));
    end
    W_eTRCA.val=[];
    for jj=1:num_class
        X0All = XTrain_fbi(YTrain==jj,:,:);
        num_of_trials=size(X0All,1);
        

        trca_X1 = squeeze(sum(X0All));
        trca_X2 = reshape(permute(X0All,[2,1,3]),num_channel, num_sampls*num_of_trials)';
        
        S=trca_X1*trca_X1'-trca_X2'*trca_X2;
        trca_X2 = trca_X2-mean(trca_X2);
        Q=trca_X2'*trca_X2;
    %     [eig_v1,eig_d1]=eig(Q\S);
        [eig_v1,eig_d1]=eig(S,Q);
        [eig_val,sort_idx]=sort(diag(eig_d1),'descend');
        eig_vec=eig_v1(:,sort_idx);
        W_eTRCA.val=[W_eTRCA.val; eig_vec(:,1)'];
    end
    yP=zeros(size(YTest));
    for trail_i=1:num_test
        XTestTemp=squeeze(XTest_fbi(trail_i,:,:));
        rowTemp = zeros(1,num_class);
        for class_i = 1:num_class
            XTestTempw=(XTestTemp'*W_eTRCA.val')';
            templateTempw=(Template{class_i}'*W_eTRCA.val')';
            rowTemp(class_i)=corr2(XTestTempw,templateTempw);
        end
        crall(trail_i,fb_i,:)=rowTemp;
    end
end
yP=zeros(size(YTest));
for i=1:num_test
    crTemp=squeeze(crall(i,:,:));
    crTemp=fb*crTemp;
    [~,yP(i)]=max(crTemp);
end
acc=sum(yP==YTest)/length(YTest);
end