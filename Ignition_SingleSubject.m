clear;
Allcell = load ('ts_menstrualcycle.mat', 'ts_menstrualcycle');

% cell array with 3 columns. Column 1 follicular, 2 pre-ovulatory, 3 luteal
% in each cell of the colum, a matrix with nregionsxtimepoints (timeseries for each subject) 
timeSeries = Allcell.ts_menstrualcycle;  

for cond = 1:3  % cond=1 follicular, cond=2 pre-ovulatory, cond=3 luteal

    % adapt parameters to YOUR data
    TR=2.25;  % Repetition Time (seconds) 225
    NSUB=60;  % total subjects
    N = 116;  % total nodes


    %%%%%%%%%%%%%

    Isubdiag = find(tril(ones(N),-1));
    nTRs = 5; % nTRs -1: TRs to compute ignition after spontanous events
    CASE = cond;
    nevents = zeros(1,N);
    FC = zeros(NSUB,N,N);

    flp = .01;              % lowpass frequency of filter
    fhi = .09;              % highpass
    delt = TR;              % sampling interval
    k = 2;                  % 2nd order butterworth filter
    fnq = 1/(2*delt);       % Nyquist frequency
    Wn = [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
    [bfilt2,afilt2] = butter(k,Wn);   % construct the filter
    ttotal = 1;

    % compute ignition for each subject
    for nsub = 1:NSUB
        nsub
        xs = timeSeries{nsub,cond};  % obtain Nregionsxtime matrix (timeseries for each subject)
        Tmax = size(xs,2); % time points
        T = 1:Tmax; 
        clear  x timeseriedata events ev1 ev2

        % obtain events for each seed
        for seed = 1:N
            x = demean(detrend(xs(seed,:)));
            timeseriedata(seed,:) = filtfilt(bfilt2,afilt2,x);    % zero phase filter the data
            Xanalytic = hilbert(demean(timeseriedata(seed,:)));
            tise = detrend(demean(timeseriedata(seed,T)));
            ev1 = tise>(std(tise)+mean(tise));
            ev2 = [0 ev1(1:end-1)];
            events(seed,:) = (ev1-ev2)>0;
        end
        SubjEvents{nsub} = events;

        %%% integration
        % obtain 'events connectivity matrix' and integration value (integ)
        % for each time point
        for t = T
            for i = 1:N
                for j = 1:N
                    phasematrix(i,j) = events(i,t)*events(j,t); 
                end
            end
            cc = phasematrix; %*Cbin;
            cc = cc-eye(N);

            [comps csize] = get_components(cc);
            integ(t) = max(csize)/N;
        end % end obtain integ

        %%%% event trigger
        nevents2 = zeros(1,N);
        % save events and integration values for nTRs after the event
        for seed = 1:N
            flag = 0;
            for t = T
                % detect first event (nevents = matrix with 1xnode and number of events in each cell)
                if events(seed,t) == 1 && flag == 0  % if events(seed,t-9) == 1 && flag == 0
                    flag = 1;
                    % events for each subject
                    nevents2(seed) = nevents2(seed)+1;
                end
                % save integration value for nTRs after the first event (nodesx(nTR-1)xevents)
                if flag > 0
                    % integration for each subject
                    IntegStim2(seed,flag,nevents2(seed)) = integ(t); 
                    flag = flag+1;
                end
                % after nTRs, set flag to 0 and wait for the next event (then, integ saved for nTRs -1 events)
                if flag == nTRs
                    flag = 0;
                end
            end
        end

        % std of the max ignition in the nTRs for each subject and for each node
        for seed = 1:N
            stdevokedinteg2(seed) = std(max(squeeze(IntegStim2(seed,:,1:nevents2(seed)))));
        end

        % std ignition across events for each subject in each node(Single Subject, S)
        stdevokedintegS(:,nsub) = stdevokedinteg2;

    end % end loop compute ignition for each subject

    save (sprintf('Ignition_CASE%d.mat', CASE), 'stdevokedintegS');
    clearvars -except condition timeSeries

end    % end loop over conditions
