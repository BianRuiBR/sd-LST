function acc = ms_ecca_fb(XTrain,YTrain,freqs,phases, XTest,YTest,num_of_signal_templates,fs,num_fb, num_harms)
%% output the accuracy of msCCA with filter bank
%   XTrain: num_train * num_channel * num_sampls * num_fb
%   YTrain: num_train * 1
%   XTest: num_test * num_channel * num_sampls * num_fb
%   YTest: num_test * 1
%   freqs: stimulus frequencies, 1*num_class
%   phases: stimulus phases, 1*num_class
%   num_of_signal_templates: number of adjacent stimuli 
%   fs: sampling frequency
%   num_fb: number of the filter bank
%   num_harms: number of harmonics

[~,id] = sort(freqs);
[~, ~, num_sampls] = size(XTrain(:,:,:,1));
class_all = unique(YTrain);
num_test = size(XTest, 1);
num_class = length(class_all);
d1=num_class;
crall3=zeros(num_test,num_fb,num_class);
fb = (1:num_fb).^(-1.25)+0.25;
for fb_i=1:num_fb
    XTrain_fbi=XTrain(:,:,:,fb_i);
    XTest_fbi=XTest(:,:,:,fb_i);

    Template=cell(1,num_class);
    ref1=cell(1,num_class);
    for class_i = 1:num_class
        XTrain_class_i = XTrain_fbi(YTrain==class_all(class_i),:,:);
        Template{class_i} = squeeze(mean(XTrain_class_i,1));
        ref1{class_i} = ref_signal_nh(freqs(class_i),fs,phases(class_i),num_sampls,num_harms);
    end
    d0=floor(num_of_signal_templates/2);
    mscca_template=cell(1,num_class);
    mscca_ref=cell(1,num_class);
    for j=1:num_class
        if j<=d0
            template_st=1;
            template_ed=num_of_signal_templates;
        elseif ((j>d0) && j<(d1-d0+1))
            template_st=j-d0;
            template_ed=j+(num_of_signal_templates-d0-1);
        else
            template_st=(d1-num_of_signal_templates+1);
            template_ed=d1;
        end
        mscca_template{j}=[];
        mscca_ref{j}=[];
        template_seq=[template_st:template_ed];

        % Concatenation of the templates (or sine-cosine references)
        for n_temp=1:num_of_signal_templates
            jj = template_seq(n_temp);
            template0=Template{id(jj)};
    %             if (is_center_std==1)
    %                 template0=template0-mean(template0,2)*ones(1,length(template0));
    %                 template0=template0./(std(template0')'*ones(1,length(template0)));
    %             end
            ref0=ref1{id(jj)};
            mscca_template{j}=[mscca_template{j};template0'];
            mscca_ref{j}=[mscca_ref{j};ref0'];
        end
    end

    Wx1=cell(1,num_class);
    Wy1=cell(1,num_class);
    for j=1:num_class
        [Wx1{j},Wy1{j},~]=canoncorr(mscca_template{j},mscca_ref{j});
    end
    for i=1:num_test
        test_signal = squeeze(XTest_fbi(i,:,:)); 
        msccaR3=zeros(1,num_class);
        for j=1:num_class    
            % ========mscca spatial filter=====
    %         [Wx1,Wy1,~]=canoncorr(mscca_template{j},mscca_ref{j});
            spatial_filter1.wx1=Wx1{j}(:,1)';
            spatial_filter1.wy1=Wy1{j}(:,1)';
            cr1=corrcoef((spatial_filter1.wx1*test_signal)',(spatial_filter1.wy1*ref1{id(j)})');
            cr2=corrcoef((spatial_filter1.wx1*test_signal)',(spatial_filter1.wx1*Template{id(j)})'); 
            msccaR3(id(j))=sign(cr1(1,2))*cr1(1,2)^2+sign(cr2(1,2))*cr2(1,2)^2; % 
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