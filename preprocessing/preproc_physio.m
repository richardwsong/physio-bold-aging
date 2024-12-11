% Script for processing physiological data from fMRI sessions
% This script processes physiological data including respiration, CO2, and cardiac measures
% from fMRI sessions, generates relevant plots, and saves processed outputs.

% Set up paths and load data
dir = '/path/to/output/';              % Output directory for processed data
D = "/path/to/physio_data/";           % Source directory containing raw physiological data
physio_QA = load('/path/to/subject_info.mat');  % Load subject information
session = 'pre';                       % Session identifier (e.g., 'pre', 'post')
id = string(physio_QA.id);            % Array of subject IDs

for i = 1:numel(id)
    try
        name = id(i);
        disp(" ");
        disp(append("starting on ", name));

        % Create subject-specific output directory
        fullpath = append(dir, name); 
        if ~isdir(fullpath)
            mkdir(fullpath);
        end

        % Load and parse physiological data JSON metadata
        fname = append(D, '/', name, '/ses-', session, '/func/', name, '_ses-', session, '_task-rest_physio.json');
        fid = fopen(fname); 
        raw = fread(fid,inf); 
        str = char(raw'); 
        fclose(fid); 
        val = jsondecode(str);
        fs_phys = 1000;                % Sampling frequency in Hz
        dt_phys = 1/fs_phys;           % Time step in seconds

        % Load physiological time series data
        fname = append(D, '/', name, '/ses-', session, '/func/');
        dat = load(append(fname, name, '_ses-', session, '_task-rest_physio.tsv'));

        % Extract scan trigger information
        tidx = find(string(val.Columns) == 'tr');
        trigger = dat(:, tidx); 
        if trigger(1) < 1
            markers = find(trigger > 1); 
        else
            markers = find(trigger < 3); 
        end
        start_scan = markers(1);           % First trigger
        end_scan = start_scan + 4200614;   % End of scan period
        TRIGGERS.start = start_scan; 
        TRIGGERS.end = end_scan; 
        TRIGGERS.tr = trigger; 

        % Extract and process respiration data
        % ------------------------------------------ %
        resp_idx = find(string(val.Columns) == 'rsp');
        resp_dat = dat(start_scan:end_scan,resp_idx);
        resp_dat = resp_dat([1:10:end]);   % Downsample by factor of 10
        resp_dat = double(resp_dat);

        % Extract and process CO2 data
        cng_idx = find(string(val.Columns) == 'cng');
        cng = dat(start_scan:end_scan,cng_idx);
        cng_ds = cng([1:10:end]);         % Downsample by factor of 10

        % Calculate time lag between CO2 and respiration using cross-covariance
        [c,lags] = xcov([cng_ds resp_dat], 24000);
        cov_vector = c(:, 2);
        index_of_zero = find(lags == 0);

        % Find optimal CO2 lag using minimum covariance
        [min_cov, min_cov_lag] = min(cov_vector(index_of_zero:end, :));
        co2_lag_ds = lags(min_cov_lag + index_of_zero);
        co2_lag = co2_lag_ds * 10;        % Scale lag back to original sampling rate

        % Re-extract CO2 with corrected lag
        cng = dat(start_scan+co2_lag:end_scan+co2_lag,cng_idx);
        cng_ds = cng([1:10:end]);

        % Find CO2 peaks with adaptive minimum peak prominence
        pks = []; 
        locs = []; 
        mpp = 10;                         % Initial minimum peak prominence
        while numel(pks) < 150 && mpp >= 0
            [pks,locs] = findpeaks(cng_ds, 'MinPeakProminence',mpp);
            mpp = mpp - 0.1; 
        end
        maxtab_c(:,1) = locs; 
        maxtab_c(:,2) = pks; 

        % Interpolate CO2 peaks to regular sampling grid
        xq = 0:2400:2400*174;            % Create regular time points
        vq1 = interp1(maxtab_c(:,1),maxtab_c(:,2),xq, 'linear', 'extrap');

        % Verify expected number of interpolated points
        disp(numel(vq1));
        assert(numel(vq1) == 175); 

        % Store respiratory data
        RESP.resp = resp_dat; 
        RESP.cng = cng; 
        RESP.cng_ds = cng_ds; 
        RESP.etco2 = vq1';
        RESP.samples = xq'; 

        % Generate and save CO2 visualization plots
        figure; 
        g1 = subplot(2, 1, 1); hold on; 
        plot(cng_ds); hold on; 
        plot(locs,pks,'r.','markersize',20);
        ylabel("cng ds")
        xlabel("physio sample #")
        title(append(name, " end-tidal CO2 peaks, lag = ", num2str(co2_lag))); 
        g2 = subplot(2, 1, 2); hold on; 
        plot(cng_ds); hold on; 
        plot(xq,vq1,'r.','markersize',20);
        ylabel("cng ds");
        xlabel("physio sample #");
        title(append(name, " end-tidal CO2 peaks interp")); 
        
        savepath = append(dir, name, '/', name, '_ses-', session, '_task-rest_physio_cng.png');
        saveas(gcf, savepath); 
    
        % Generate and save respiration/CO2 comparison plot
        figure; 
        plot((resp_dat-mean(resp_dat))/mean(resp_dat)); 
        hold on; 
        plot((cng_ds-mean(cng_ds))/mean(cng_ds)/15 + .2)
        savepath = append(dir, name, '/', name, '_ses-', session, '_task-rest_physio_cng_resp.png');
        saveas(gcf, savepath); 

        % Process Cardiac (PPG) Data
        % ------------------------------------------ %
        ppg_idx = find(string(val.Columns) == 'ppg');
        ppg = dat(start_scan:end_scan,ppg_idx);
        card_dat = ppg([1:10:end]);       % Downsample by factor of 10
        
        % Apply bandpass filter for cardiac signal
        fcut_BPF = [0.5,2];               % Cutoff frequencies in Hz
        Fn = fs_phys/2;                   % Nyquist frequency
        Wn = fcut_BPF/Fn; 
        Nb = 2;                           % Filter order
        [B, A] = butter(Nb,Wn);
        card_bpf =  filtfilt(B,A,double(card_dat));

        % Detect cardiac peaks
        card_rng = iqr(card_bpf); 
        minHeight = 0.05*card_rng;
        [pks,locs] = findpeaks(card_bpf,'minpeakheight',minHeight);
        maxtab_c(:,1) = locs; 
        maxtab_c(:,2) = pks;

        % Calculate cardiac timing and heart rate metrics
        card_trig_samples = locs;
        card_trig_times = card_trig_samples*dt_phys;
        IBI = diff(card_trig_times);      % Inter-beat intervals
        HR = (1./diff(card_trig_times))*60; % Heart rate in BPM

        % Generate and save cardiac data visualization
        figure; 
        g1 = subplot(3,1,1); hold on;
        plot(card_bpf); 
        plot(maxtab_c(:,1),maxtab_c(:,2),'r.','markersize',20);
        legend('cardiac data (band-pass filt)','detected beats');
        xlabel('physio sample #');
        
        g2 = subplot(3,1,2); hold on;
        plot(card_trig_samples(1:end-1),HR,'b'); 
        ylabel('beats per min'); 
        title('cardiac rate');
        hold on;
        plot(card_trig_samples(1:end-1),HR,'c.');
        xlabel('physio sample #');
        linkaxes([g1,g2],'x');
        
        g3 = subplot(3,1,3); hold on;
        plot(IBI); 
        xlabel('index'); 
        ylabel('IBI');
        drawnow;

        save_path_2 = append(dir, name, '/', name, '_ses-', session, '_task-rest_physio_cardiac.png');
        saveas(gcf, save_path_2); 
        
        % Create output structure with processed cardiac data
        OUT_p.IBI_raw = IBI;
        OUT_p.HR_raw = HR;
        OUT_p.card_trig_times_s = card_trig_times;
        OUT_p.card_trig_samples = card_trig_samples;
        OUT_p.card_dat = card_dat;
        OUT_p.card_bpf = card_bpf;
        OUT_p.card_pks = maxtab_c;
        OUT_p.dt_phys = dt_phys;
        OUT_p.cng = cng_ds;
        OUT_p.resp = resp_dat; 
        
        % Remove outlier beats and generate clean heart rate metrics
        [OUT_p.IBI_clean,OUT_p.HR_clean,outliers_IBI, outliers_HR] = despike_hr(OUT_p);
        OUT_p.outlier_IBI = outliers_IBI;
        OUT_p.outlier_HR = outliers_HR; 
        if ~isempty(outliers_IBI) || ~isempty(outliers_HR)
            save_path_4 = append(dir, name, '/', name, '_ses-', session, '_task-rest_physio_outliers.png');
            saveas(gcf, save_path_4);       
        end
        
        % Generate fMRI physiological regressors
        REGS = build_fmri_regs(OUT_p,2.4,175);
        save_path_3 = append(dir, name, '/', name, '_ses-', session, '_task-rest_physio_regs.png');
        saveas(gcf, save_path_3); 

        % Save all processed data
        savePath = append(dir, name, '/', name, '_ses-', session, '_task-rest_physio','_physOUT.mat');
        save(savePath,'OUT_p','RESP','REGS', 'TRIGGERS');

        close all; 
    catch e
        % Error handling: save trigger information if processing fails
        disp(append('failed on ', name)); 
        disp(e.message);
        TRIGGERS.start = start_scan; 
        TRIGGERS.end = end_scan; 
        TRIGGERS.tr = trigger; 
        savePath = append(dir, name, '/', name, '_ses-', session, '_task-rest_physio','_physOUT.mat');
        save(savePath,'TRIGGERS');
    end
end

function [IBI_clean,HR_clean,outliers_IBI, outliers_HR] = despike_hr(IN_p)
% Removes outliers from heart rate and IBI time series
%
% Parameters:
%   IN_p: Structure containing raw cardiac data
%
% Returns:
%   IBI_clean: Cleaned inter-beat interval time series
%   HR_clean: Cleaned heart rate time series
%   outliers_IBI: Indices of IBI outliers
%   outliers_HR: Indices of HR outliers

    outliers_IBI = [find(IN_p.IBI_raw >= mean(IN_p.IBI_raw) + 2.5*std(IN_p.IBI_raw)); 
                    find(IN_p.IBI_raw <= mean(IN_p.IBI_raw) - 2.5*std(IN_p.IBI_raw))]; 
    outliers_HR = [find(IN_p.HR_raw >= mean(IN_p.HR_raw) + 2.5*std(IN_p.HR_raw)); 
                    find(IN_p.HR_raw <= mean(IN_p.HR_raw) - 2.5*std(IN_p.HR_raw))]; 
    
    figure;
    subplot(211);
    IBI_clean = interp_ts(IN_p.IBI_raw, outliers_IBI, 1);
    subplot(212);
    HR_clean = interp_ts(IN_p.HR_raw, outliers_HR, 1);
end

function REGS = build_fmri_regs(IN_p,TR,nframes)
% Generates physiological regressors for fMRI analysis
%
% Parameters:
%   IN_p: Structure containing processed physiological data
%   TR: Repetition time of fMRI acquisition in seconds
%   nframes: Number of fMRI frames
%
% Returns:
%   REGS: Structure containing physiological regressors aligned to fMRI timing
    
    dt_phys = IN_p.dt_phys;
    resp = IN_p.resp;
    card_dat = IN_p.card_dat;
    card_bpf = IN_p.card_bpf;
    IBI_clean = IN_p.IBI_clean;
    
    % Calculate timing for IBI values
    t_ibi = 0.5*(IN_p.card_trig_times_s(2:end) + IN_p.card_trig_times_s(1:end-1));
    assert(length(t_ibi)==length(IBI_clean))
    
    % Set up fMRI timing parameters
    Twin = 6;                          % Window size in seconds
    TR_s = TR;
    t_fmri = (TR_s/2)+[0:TR_s:TR_s*nframes-TR_s];
    
    % Initialize regressor arrays
    rv = [];                           % Respiratory variation
    pa = [];                           % Pulse amplitude
    hr = [];                           % Heart rate
    pa_bpf = [];                       % Filtered pulse amplitude
    
    % Calculate physiological metrics for each fMRI frame
    for kk=1:nframes
        t = t_fmri(kk);
        
        % Calculate heart rate metrics
        t1 = max(0,t-Twin*0.5);
        t2 = min(TR_s*nframes,t+Twin*0.5);
        inds = intersect(find(t_ibi<=t2),find(t_ibi>=t1));
        hr(kk) = (60./median(IBI_clean(inds)));
        hrv(kk) = sqrt(mean(diff(IBI_clean(inds)).^2)); % RMSSD
        
        % pulse amplitude
        % ---------------------- %
        if length(resp)~=length(card_dat)
            error('resp & card sampled at different rates');
        else
            np = length(resp);
        end
        
        % window (in samples)
        i1 = max(1,floor((t - Twin*0.5)/dt_phys)); 
        i2 = min(np, floor((t + Twin*0.5)/dt_phys));
        pa(kk) = std(card_dat(i1:i2));
        pa_bpf(kk) = std(card_bpf(i1:i2));

        % respiration variation
        % ---------------------- %
        rv(kk) = std(resp(i1:i2));
    end
    
    % regressors for fmri
    % ---------------------- %
    REGS.rv = rv(:);
    REGS.pa = pa(:);
    REGS.pa_bpf = pa_bpf(:);
    REGS.hr = hr(:);
    REGS.hrv = hrv(:);
    
    figure;
    subplot(411); 
    plot(hr); title('heart rate');
    subplot(412);
    plot(hrv); title('heart rate variability (rmssd)');
    subplot(413);
    plot(pa); title('pulseox amplitude');
    hold on;
    plot(pa_bpf,'r'); title('pulseox amplitude - cardbpf');
    subplot(414);
    plot(rv); title('respiratory variation');
    xlabel('fmri TR');
end