function acc = ttCCA(XSource, YSource, XTest,YTest,freqs,fs,num_fb, num_harms)
%% output the accuracy of ttCCA
%  XSource: cell(1，num_Source), XSource{1}: num_trial * num_channel * num_sampls * num_fb
%   YSource: cell(1，num_Source), YSource{1}: num_trial * 1
%   XTest: num_test * num_channel * num_sampls * num_fb
%   YTest: num_test * 1
%   freqs: stimulus frequencies, 1*num_class
%   fs: sampling frequency
%   num_fb: number of the filter bank
%   num_harms: number of harmonics

num_Source = length(XSource);
[~, num_channel, num_sampls] = size(XTest(:,:,:,1));
class_all = unique(YTest);
num_class = length(class_all);
num_test = size(XTest(:,:,:,1), 1);
Yref=cell(1,num_class);
for class_i = 1:num_class
    Yref{class_i} = ref_signal_nh(freqs(class_i),fs,0,num_sampls,num_harms);
end
fb = (1:num_fb).^(-1.25)+0.25;
crall = zeros(num_test,num_fb,num_class);
for fb_i =1:num_fb
    XTest_fbi=XTest(:,:,:,fb_i);
    TamplateSource = cell(1,num_class);
    for class_i =1:num_class
        TamplateSource{class_i} = zeros(num_channel, num_sampls);
        for sub_i =1:num_Source
            TamplateSource{class_i} = TamplateSource{class_i} + squeeze(mean(XSource{sub_i}(YSource{sub_i} == class_i,:,:,fb_i),1));
        end
        TamplateSource{class_i} = TamplateSource{class_i} / num_Source;
    end
    for test_i = 1:num_test
        test_signal = squeeze(XTest_fbi(test_i,:,:)); 
        for j=1:num_class
            [Wx,~,row2] = canoncorr(test_signal',Yref{j}');
            row2 = row2(1);
            [Wxt,~,~] = canoncorr(TamplateSource{j}',Yref{j}');
            row1 = corrcoef((Wx(:,1)'*test_signal)',(Wx(:,1)'*TamplateSource{j})');
            row1 = row1(1,2);
            row3 = corrcoef((Wxt(:,1)'*test_signal)',(Wxt(:,1)'*TamplateSource{j})');
            row3 = row3(1,2);
            crall(test_i,fb_i,j)=row1 + row2 + row3;
        end
    end
end

yP=zeros(size(YTest));
for i=1:num_test
    crTemp=squeeze(crall(i,:,:));
    crTemp=fb*crTemp;
    [~,yP(i)]=max(crTemp);
end
acc=sum(yP==YTest)/num_test;
end

function refSignal=ref_signal_nh(f,fs,phase,tlen,num_of_harmonics)
% constract first and second harmonic wave
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