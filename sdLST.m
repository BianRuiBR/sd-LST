function acc = sdLST(XTarget,YTarget,XSource,YSource,XTest,YTest,freqs,phases,fs,num_fb, num_harms)
%% output the accuracy of sdLST with filter bank
%   XTarget: num_target * num_channel * num_sampls * num_fb
%   YTarget: num_target * 1
%   XSource: cell(1，num_Source), XSource{1}: num_trial * num_channel * num_sampls * num_fb
%   YSource: cell(1，num_Source), YSource{1}: num_trial * 1
%   XTest: num_test * num_channel * num_sampls * num_fb
%   YTest: num_test * 1
%   freqs: stimulus frequencies, 1*num_class
%   phases: stimulus phases, 1*num_class
%   fs: sampling frequency
%   num_fb: number of the filter bank
%   num_harms: number of harmonics

num_Source = length(XSource);
source_Num = size(XSource{1},1);
[~, num_channel, num_sampls] = size(XTest(:,:,:,1));
targetClass = unique(YTarget);
class_all = unique(YSource{1});
num_test = size(XTest(:,:,:,1), 1);
num_class = length(class_all);
Yref=cell(1,num_class);
for class_i = 1:num_class
    Yref{class_i} = ref_signal_nh(freqs(class_i),fs,phases(class_i),num_sampls,num_harms);
end
YSourceAll = zeros(num_Source * source_Num, 1);
crall=zeros(num_test,num_fb,num_class);
for fb_i=1:num_fb
    XSourceAll = zeros(num_Source * source_Num, num_channel, num_sampls);
    XTest_fbi=XTest(:,:,:,fb_i);
    XTarget_fbi=XTarget(:,:,:,fb_i);
    XTarget_fbi_tamp = [];
    for tar_i = 1:length(targetClass)
        XTarget_fbi_tamp = [XTarget_fbi_tamp,squeeze(mean(XTarget_fbi(YTarget == targetClass(tar_i),:,:), 1))];
    end
    for sub_i = 1:num_Source
        XSource_fbi_tamp = [];
        for tar_i = 1:length(targetClass)
            XSource_fbi_tamp = [XSource_fbi_tamp,squeeze(mean(XSource{sub_i}(YSource{sub_i} == targetClass(tar_i),:,:,fb_i),1))];
        end
        P = XTarget_fbi_tamp*XSource_fbi_tamp'/(XSource_fbi_tamp*XSource_fbi_tamp');
        for trial_i = 1:size(XSource{sub_i},1)
            tmp = squeeze(XSource{sub_i}(trial_i,:,:,fb_i));
            tmp = P*tmp;
            XSourceAll((sub_i-1)*source_Num + trial_i, :, :) = tmp;
            YSourceAll((sub_i-1)*source_Num + trial_i) = YSource{sub_i}(trial_i);
        end
    end
    XSource_fbi_tamp = [];
    for tar_i = 1:length(targetClass)
        XSource_fbi_tamp = [XSource_fbi_tamp,squeeze(mean(XSourceAll(YSourceAll == targetClass(tar_i),:,:),1))];
    end
    P = XTarget_fbi_tamp*XSource_fbi_tamp'/(XSource_fbi_tamp*XSource_fbi_tamp');
    for trial_i =1:length(YSourceAll)
        tmp = squeeze(XSourceAll(trial_i,:,:));
        tmp = P*tmp;
        XSourceAll(trial_i, :, :) = tmp;
    end
    Tamplate = cell(1,num_class);
    for class_i = 1:num_class
        Tamplate{class_i} = squeeze(mean(XSourceAll(YSourceAll == class_i, :, :)));
    end
    W_eTRCA = [];
    for jj=1:num_class
        X0All = XSourceAll(YSourceAll==jj,:,:);
        num_of_trials=size(X0All,1);
        trca_X1 = squeeze(sum(X0All));
        trca_X2 = reshape(permute(X0All,[2,1,3]),num_channel, num_sampls*num_of_trials)';
        
        S=trca_X1*trca_X1'-trca_X2'*trca_X2;
        trca_X2 = trca_X2-mean(trca_X2);
        Q=trca_X2'*trca_X2;
        
        [eig_v1,eig_d1]=eig(S,Q);
        [~,sort_idx]=sort(diag(eig_d1),'descend');
        eig_vec=eig_v1(:,sort_idx);
        W_eTRCA=[W_eTRCA; eig_vec(:,1)'];
    end
    w_ecca = cell(1, num_class);
    for class_i = 1:num_class
        [w_ecca{class_i},~,~]=canoncorr(Tamplate{class_i}',Yref{class_i}');
    end
    for test_i = 1:num_test
        XTestTemp=squeeze(XTest_fbi(test_i,:,:));
        for class_i = 1:num_class
            [Wx1,Wy1,~]=canoncorr(XTestTemp',Yref{class_i}');
            row1 = corrcoef((Wx1(:,1)'*XTestTemp)',(Wy1(:,1)'*Yref{class_i})');
            row1 = row1(1,2);
            row2 = corrcoef((Wx1(:,1)'*XTestTemp)',(Wx1(:,1)'*Tamplate{class_i})');
            row2 = row2(1,2);
            row3 = corrcoef((w_ecca{class_i}(:,1)'*XTestTemp)',(w_ecca{class_i}(:,1)'*Tamplate{class_i})');
            row3 = row3(1,2);
            row0 = corr2((XTestTemp'*W_eTRCA')',(Tamplate{class_i}'*W_eTRCA')');
            crall(test_i,fb_i,class_i) = sign(row0) * row0^2 + (sign(row1) * row1^2 + sign(row2) * row2^2 + sign(row3) * row3^2) / 3;
        end
    end
end
fb = (1:num_fb).^(-1.25)+0.25;
yP=zeros(size(YTest));
for i=1:num_test
    crTemp=squeeze(crall(i,:,:));
    crTemp=fb*crTemp;
    [~,yP(i)]=max(crTemp);
end
acc=sum(yP==YTest)/length(YTest);
end