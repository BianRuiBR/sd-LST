function acc = stCCA(XTarget,YTarget,XSource,YSource,XTest,YTest,freqs,phases,fs, num_fb,num_harms)
%% output the accuracy of stCCA with filter bank
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

warning off;
num_Source = length(XSource);
[~, ~, num_sampls] = size(XTest(:,:,:,1));

targetClass = unique(YTarget);
K = length(targetClass);
class_all = unique(YSource{1});
num_test = size(XTest(:,:,:,1), 1);
num_class = length(class_all);
Yref=cell(1,num_class);
for class_i = 1:num_class
    Yref{class_i} = ref_signal_nh(freqs(class_i),fs,phases(class_i),num_sampls,num_harms);
end
fb = (1:num_fb).^(-1.25)+0.25;
crall3=zeros(num_test,num_fb,num_class);
for fb_i=1:num_fb
    XTarget_fbi=XTarget(:,:,:,fb_i);
    XTest_fbi=XTest(:,:,:,fb_i);
    XTarget_fbi_tamp = [];
    XtargetAll = [];
    YrefTargetAll = [];
    for class_i =1:K
        tmp = squeeze(mean(XTarget_fbi(YTarget == targetClass(class_i),:,:), 1));
        XTarget_fbi_tamp(class_i,:,:) = tmp;
        XtargetAll = [XtargetAll,tmp];
        YrefTargetAll = [YrefTargetAll, Yref{targetClass(class_i)}];
    end
    [Wx1target,Wy1target,~]=canoncorr(XtargetAll',YrefTargetAll');
    Wx1target = Wx1target(:,1);
    Wy1target = Wy1target(:,1);
    bbb = [];
    for class_i =1:K
        bbb = [bbb;(Wx1target'*squeeze(XTarget_fbi_tamp(class_i,:,:)))'];
    end
    tamplateSourceAll = cell(num_Source,num_class);
    Wx1Source = cell(1,num_Source);
    AAA = [];
    for sub_i = 1:num_Source
        aaa = [];
        for class_i =1:num_class
            tmp = squeeze(mean(XSource{sub_i}(YSource{sub_i} == class_i,:,:,fb_i),1));
            tamplateSourceAll{sub_i, class_i} = tmp;
        end
        XSourceAll = [];
        for class_i =1:K
            XSourceAll = [XSourceAll,tamplateSourceAll{sub_i, targetClass(class_i)}];
        end
        [Wx1Source{sub_i},~,~]=canoncorr(XSourceAll',YrefTargetAll');
        Wx1Source{sub_i} = Wx1Source{sub_i}(:,1);
        for class_i =1:K
            aaa = [aaa;(Wx1Source{sub_i}'*tamplateSourceAll{sub_i, targetClass(class_i)})'];
        end
        AAA = [AAA,aaa];
    end
    www = (AAA'*AAA)\AAA'*bbb;
    tamplateForTarget = cell(1,num_class);
    for class_i = 1:num_class
        tamplateForTarget{class_i} = zeros(num_sampls,1);
        for sub_i = 1: num_Source
            tamplateForTarget{class_i} =  tamplateForTarget{class_i} + www(sub_i)*(Wx1Source{sub_i}'*tamplateSourceAll{sub_i, class_i})';
        end
    end
    for i=1:num_test
        test_signal = squeeze(XTest_fbi(i,:,:)); 
        msccaR3=zeros(1,num_class);
        for j=1:num_class    
            cr2=corrcoef((Wx1target'*test_signal)',tamplateForTarget{j});
            cr1=corrcoef((Wx1target'*test_signal)',(Wy1target'*Yref{j})'); 
            msccaR3(j)=sign(cr1(1,2))*cr1(1,2)^2+sign(cr2(1,2))*cr2(1,2)^2; % 
        end
        crall3(i,fb_i,:)=msccaR3;
    end
end

yP3=zeros(size(YTest));
for i=1:num_test
    crTemp=squeeze(crall3(i,:,:));
    crTemp=fb*crTemp;
    [~,yP3(i)]=max(crTemp);
end
acc=sum(yP3==YTest)/num_test;

end

function refSignal=ref_signal_nh(f,fs,phase,tlen,num_of_harmonics)
%% constract first and second harmonic wave
p=tlen;
% fs=250;
% TP=1/fs:1/fs:p/fs;
TP=0:1/fs:p/fs-1/fs;
refSignal=[];
for h=1:num_of_harmonics
    Sinh1=sin(2*pi*h*f*TP+h*phase);
    Cosh1=cos(2*pi*h*f*TP+h*phase);
%     Sinh2=sin(2*pi*2*f*TP+2*phase);
%     Cosh2=cos(2*pi*2*f*TP+2*phase);
%     Sinh3=sin(2*pi*3*f*TP+3*phase);
%     Cosh3=cos(2*pi*3*f*TP+3*phase);
    refSignal=[refSignal;Sinh1;Cosh1;];
end
end