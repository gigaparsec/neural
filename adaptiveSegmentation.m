%% Author: Jack Tchimino

% Adaptive segmentation function. This function will divide the signal into
% segments, using the spectral error measure criterion and the ACF criterion. 
% follows the algorithm in Rangayyan's Biomedical Signal Analysis, page 489

%%
% threshold: the SEM threshold at which the signal will be segmented. P is
% the order of the AR model that willbe used. plotYN is a "Y" or "N" input,
% which makes the function generate the plots or not. method can be ACF or
% SEM. N: number of samples in window

function segments = adaptiveSegmentation(data,threshold,N,P,plotYN,method)

% the n counter is different than the center of the window n0!!!!

    signal = data.signal;
    fs = data.SamplingFrequency;
    tV = data.TimeVector;

    n0 = N+P+1;
    axSt=n0;
    NS = length(signal);
%     phi_e = zeros(NS,2*P+1);
    phi_e = zeros(NS,P);
%     
if strcmp(method,'SEM')==1
    init_n=1;
    seg_num=0;
    while n0<NS-N-P %<------------------------- loop start
        % define initial window and find its ACF
        initWindow = signal(n0-N:n0+N);
        signalACF = xcorr(initWindow,P);
        % ar model
        temp = iddata(double(initWindow),[],1/fs);
        ARmod = ar(temp,P);
        % find prediction error
        windForPE = signal(n0-N-P:n0+N);
        predError = pe(ARmod,(windForPE),N);
        % compute short time ACF phi_e(n,m) for PE
        m0=P+1;
        for m=0:1:P
            factor = 1/(2*N + 1);
            temp = 0;
            disp(['-N=', num2str(-N),', N-m=',num2str(N-m)])
    %         pause
            for k=-N:1:N-m %not N-m
                disp(['n=',num2str(n0),' m=',num2str(m),', k=',num2str(k),', n+k=',num2str(n0+k),', n+k+m=',num2str(n0+k+m)])
%                 temp = temp+ predError(n0+k)*predError(n0+k+m);
                temp = temp+ predError(axSt+k)*predError(axSt+k+m);
    %             n+k
    %             n+k+m
            end
            phi_e(n0,m+1) = temp*factor;
        end
        
        fixedWindowPhi_e = phi_e(n0,:);
        
        SEM=0;
        
%         while SEM<threshold && n<NS-N-P %<----------- loop2 start
        n0=n0+1;
        movingWindowPhi_e = zeros(size(fixedWindowPhi_e));
        for m=1:1:P+1
%             movingWindowPhi_e(n0,m) = phi_e(n0-1,m)+predError(n0+N)*predError(n0+N-m)-predError(n0-N-1)*predError(n0-N-1-m);
            a = phi_e(n0-1,m);
            b = predError(axSt-1+N);
            c = predError(axSt+N-1-m);
            d = predError(axSt-N);
            e = predError(1+axSt-N-m);
            phi_e(n0,m) = a+b*c-d*e;
        end
        
        fac1 = ((fixedWindowPhi_e(1)/phi_e(n0,1)) -1)^2;
%         fac2 = 
        summed=0;
        for k=1:1:P
            summed=summed+(phi_e(n0,k)/phi_e(n0,1))^2;
        end
        fac2 = 2*summed;
        SEM=fac1+fac2;
        SpectralErrorMeasure(n0) = SEM;
        
        if SEM>threshold
                flag=1;
            elseif SEM<threshold && n0>NS-N-P
                flag=2;
        end
        
        
%         end %<--------------------------------------- loop2 end
        
        if flag==2
            disp('No further segmentation possible')
            segments = cell(0);
            return;
        elseif flag==1
            disp(['Segment found at Sample n=',num2str(n0)]);
            seg_num=seg_num+1;
            segments{seg_num}=[init_n n0-1];
            initWindow = signal(n0-N:n0+N); % the reference window is redefined as the last sliding window
            n0=n0+1;
            flag=0;
        end

        
%         segments = cell(0);
    end %<------------------------------------- loop end


    %step1
%     
%     initWindow = signal(n-N:n+N);
%     
%     while n<NS-N-P
%         signalACF = xcorr(initWindow,P);
%         %step2
%         temp = iddata(double(initWindow),[],1/fs);
%         ARmod = ar(temp,P);
%         %step3
%         wind = signal(n-N-P:n+N+P);% not n+N
%         predError = pe(ARmod,(wind),N);
%     %     predErrACF(n,:) = xcorr(predError,P); % use instead of eq.8.22
% 
%         % eq 8.22
%         m0=P+1;
%         for m=-P:1:P
%             factor = 1/(2*N + 1);
%             temp = 0;
%             disp(['-N=', num2str(-N),', N-m=',num2str(N-m)])
%     %         pause
%             for k=-N:1:N %not N-m
%                 disp(['n=',num2str(n),' m=',num2str(m),', k=',num2str(k),', n+k=',num2str(n+k),', n+k+m=',num2str(n+k+m)])
%                 temp = temp+ predError(n+k)*predError(n+k+m);
%     %             n+k
%     %             n+k+m
%             end
%             phi_e(n,m0+m) = temp*factor;
%         end
% 
%         %step4
%         fixedWindowPhi(n+1,:) = phi_e(n,:);
% 
%         % begin loop
%         SEM=0;
%         n=n+1;
%         while SEM<threshold && n<NS-N-P
%             movingWindowPhi = zeros(size(fixedWindowPhi));
%             for m=0:1:P
%                 movingWindowPhi(n,m0+m) = phi_e(n-1,m0+m)+predError(n+N)*predError(n+N-m)-predError(n-N-1)*predError(n-N-1-m);
%             end
%             %SEM
%             fac1 = ((fixedWindowPhi(n,m0)/movingWindowPhi(n,m0)) -1)^2;
%     %         fac2 = 
%             summed=0;
%             for k=-P:1:P
%                 summed=summed+(movingWindowPhi(n,m0+k)/movingWindowPhi(n,m0))^2;
%             end
%             fac2 = 2*summed;
%             SEM=fac1+fac2;
%             SpectralErrorMeasure(n) = SEM;
% 
%             n = n+1;
% 
%             if SEM>threshold
%                 flag=1;
%             elseif SEM<threshold && n>NS-N-P
%                 flag=2;
%             end
%         end
%         if flag==2
%             disp('No further segmentation possible')
%             segments = cell(0);
%             return;
%         elseif flag==1
%             disp(['Segment found at Sample n=',num2str(n)]);
%             seg_num=seg_num+1;
%             segments{seg_num}=[init_n n];
%             initWindow = signal(n-N:n+N); % the reference window is redefined as the last sliding window
%             n=n+1;
%             flag=0;
%         end
%     end
%     
    if isempty(segments)==false && strcmp(plotYN,'Y')==1
        figure
        subplot(2,1,1)
        plot(tv,signal)
        ax=gca;
        hold all
        for i=1:length(segments)
            seg_time = segments{i}(1);
            h1 = line([seg_time seg_time],ax.YLim);
            h1.Color='r';
            uistack(h1);
        end
        ylabel('Voltage');
        subplot(2,1,2)
        plot(tV,SEM)
        ax=gca;
        hold all
        for i=1:length(segments)
            seg_time = segments{i}(1);
            h1 = line([seg_time seg_time],ax.YLim);
            h1.Color='r';
            uistack(h1);
        end
        ylabel('Spectral Error Measure');xlabel('Time [s]');
        suptitle('Segmentation with Spectral Error Measure')
    end
elseif strcmp(method,'ACF')==1
        
else
     disp("Error: Try again");
    return;
end   