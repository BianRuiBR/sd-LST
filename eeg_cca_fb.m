function acc = eeg_cca_fb(eeg, yTest, list_freqs, fs, num_harms, num_fb)
%% output the accuracy of sCCA with filter bank， i.e. fbCCA
%   eeg: num_test * num_channel * num_sampls * num_fb
%   yTest: num_test * 1
%   list_freqs: stimulus frequencies, 1*num_class
%   fs: sampling frequency
%   num_harms: 
%   num_fb: number of the filter bank
%   num_harms: number of harmonics

[num_targs, ~, num_smpls] = size(eeg(:,:,:,1));
y_ref = cca_reference(list_freqs, fs, num_smpls, num_harms);
num_class = length(list_freqs);
r=zeros(1,num_class);
results=zeros(num_targs,1);
fb = (1:num_fb).^(-1.25)+0.25;
crall=zeros(num_targs,num_fb,num_class);
for fb_i=1:num_fb
    XTest_fbi=eeg(:,:,:,fb_i);
    for targ_i = 1:1:num_targs
        temp = XTest_fbi(targ_i, :, :);
        sizetemp=size(temp);
        test_tmp = reshape(temp,sizetemp(2:3));
        for class_i = 1:1:num_class
            refdata = squeeze(y_ref(class_i, :, :));
            [~,~,r_tmp] = canoncorr(test_tmp', refdata');
            crall(targ_i,fb_i,class_i) = r_tmp(1,1);
        end % class_i
    end
end % targ_i
yP=zeros(size(yTest));
for i=1:num_targs
    crTemp=squeeze(crall(i,:,:));
    crTemp=fb*crTemp;
    [~,yP(i)]=max(crTemp);
end
acc=sum(yP==yTest)/length(yTest);
end

function [ y_ref ] = cca_reference(list_freqs, fs, num_smpls, num_harms)
if nargin < 3 
    error('stats:cca_reference:LackOfInput',...
        'Not enough input arguments.');
end

if ~exist('num_harms', 'var') || isempty(num_harms), num_harms = 3; end

num_freqs = length(list_freqs);
tidx = (1:num_smpls)/fs;
for freq_i = 1:1:num_freqs
    tmp = [];
    for harm_i = 1:1:num_harms
        stim_freq = list_freqs(freq_i);
        tmp = [tmp;...
            sin(2*pi*tidx*harm_i*stim_freq);...
            cos(2*pi*tidx*harm_i*stim_freq)];
    end % harm_i
    y_ref(freq_i, 1:2*num_harms, 1:num_smpls) = tmp;
end % freq_i
end