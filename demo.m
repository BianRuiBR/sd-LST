%% a test demo for sdLST, data have been preprocessed, only one source subject and one target subject, so the results are not representative. 
rng(0);
fs = 256;
sub_num = 2;
freq_num = 12;
num_block = 15;
num_trial=180;
num_fb=5;
num_harms = 5;
freq_phase = load('data/Freq_Phase.mat');
list_freqs = freq_phase.freqs;
[list_freqs_sort, list_freqs_sort_id] = sort(list_freqs);
list_phases = freq_phase.phases;
list_phases_sort = list_phases(list_freqs_sort_id);

data = load(['data/datademo.mat']);
XAll = data.XAll;
YAll = data.YAll;

acc_sdLST = zeros(freq_num, sub_num,num_block);
acc_stCCA = zeros(freq_num, sub_num,num_block);

for freq_i = 6:freq_num
    K = 1:freq_i;
    K = floor(12*(2*K - 1)/(2*freq_i)) + 1;
    parfor data_i =1:sub_num
        XTargetAll = XAll{data_i};
        YTargetAll = YAll{data_i};
        XSource = XAll;
        XSource(data_i)=[];
        YSource= YAll;
        YSource(data_i)=[];
        for block_i =1:num_block
            targetId = block_i*12-11:block_i*12;
            testId = 1:num_trial;
            testId(targetId)=[];
            XTarget = XTargetAll(targetId,:,:,:);
            YTarget = YTargetAll(targetId);
            XTest = XTargetAll(testId,:,:,:);
            YTest = YTargetAll(testId);
            XTarget1 = XTarget(K,:,:,:);
            YTarget1 = YTarget(K);
            acc_sdLST(freq_i,data_i,block_i)=sdLST(XTarget1,YTarget1,XSource,YSource,XTest,YTest,list_freqs_sort,list_phases_sort,fs,num_fb,num_harms);
            acc_stCCA(freq_i,data_i,block_i)=stCCA(XTarget1,YTarget1,XSource,YSource,XTest,YTest,list_freqs_sort,list_phases_sort,fs, num_fb,num_harms);
        end
    end
end
