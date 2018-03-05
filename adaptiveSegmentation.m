%% Author: Jack Tchimino

% Adaptive segmentation function. This function will divide the signal into
% segments, using the spectral error measure criterion and the ACF criterion. 
% follows the algorithm in Ragnayyan's Biomedical Signal Analysis, page 489

%%
% threshold: the SEM threshold at which the signal will be segmented. P is
% the order of the AR model that willbe used. plotYN is a "Y" or "N" input,
% which makes the function generate the plots or not. method can be ACF or
% SEM. N: number of samples in window

function segments = adaptiveSegmentation(data,threshold,N,P,plotYN,method)

    signal = data.signal;
    fs = data.Samplingfrequency;
    tV = data.TimeVector;

    n = N+P+1;
    NS = length(signal);
    phi_e = zeros(NS,2*P+1);
    
if strcmp(method,'SEM')==1
    init_n=1;
    seg_num=0;
%     while n<NS-N-P
    %step1
    
    initWindow = signal(n-N:n+N);
    
    while n<NS-N-P
        signalACF = xcorr(initWindow,P);
        %step2
        temp = iddata(double(initWindow),[],1/fs);
        ARmod = ar(temp,P);
        %step3
        wind = signal(n-N-P:n+N+P);% not n+N
        predError = pe(ARmod,(wind),N);
    %     predErrACF(n,:) = xcorr(predError,P); % use instead of eq.8.22

        % eq 8.22
        m0=P+1;
        for m=-P:1:P
            factor = 1/(2*N + 1);
            temp = 0;
            disp(['-N=', num2str(-N),', N-m=',num2str(N-m)])
    %         pause
            for k=-N:1:N %not N-m
                disp(['n=',num2str(n),' m=',num2str(m),', k=',num2str(k),', n+k=',num2str(n+k),', n+k+m=',num2str(n+k+m)])
                temp = temp+ predError(n+k)*predError(n+k+m);
    %             n+k
    %             n+k+m
            end
            phi_e(n,m0+m) = temp*factor;
        end

        %step4
        fixedWindowPhi = phi_e(n,:);

        % begin loop
        SEM=0;
        n=n+1;
        while SEM<threshold && n<NS-N-P
            movingWindowPhi = zeros(size(fixedWindowPhi));
            for m=-P:1:P
                movingWindowPhi(n,m0+m) = phi_e(n-1,m0+m)+predError(n+N)*predError(n+N-m)-predError(n-N-1)*predError(n_N-1-m);
            end
            %SEM
            fac1 = ((fixedWindowPhi(n,m0)/movingWindowPhi(n,m0)) -1)^2;
    %         fac2 = 
            summed=0;
            for k=-P:1:P
                summed=summed+(movingWindowPhi(n,m0+k)/movingWindowPhi(n,m0))^2;
            end
            fac2 = 2*summed;
            SEM=fac1+fac2;
            SpectralErrorMeasure(n) = SEM;

            n = n+1;

            if SEM>threshold
                flag=1;
            elseif SEM<threshold && n>NS-N-P
                flag=2;
            end
        end
        if flag==2
            disp('No further segmentation possible')
            segments = cell(0);
            return;
        elseif flag==1
            disp(['Segment found at Sample n=',num2str(n)]);
            seg_num=seg_num+1;
            segments{seg_num}=[init_n n];
            initWindow = signal(n-N:n+N);
            n=n+1;
            flag=0;
        end
    end
    
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