emgdata = LoadEmgDataFile('SubjectData/Subject32/MVC.txt');

% raw data
emg_br_raw = emgdata(:,1);
emg_bi_raw = emgdata(:,2);
emg_tlong_raw = emgdata(:,3);
emg_pd_raw = emgdata(:,5);
emg_ad_raw = emgdata(:,6);

% rectified data
emg_br_rec = abs(emg_br_raw);
emg_bi_rec = abs(emg_bi_raw);
emg_tlong_rec = abs(emg_tlong_raw);
emg_pd_rec = abs(emg_pd_raw);
emg_ad_rec = abs(emg_ad_raw);

% 2nd order butterworth filter
samplerate = 2000; %hz
cutoff = 4;
[b,a] = butter(2, cutoff/(samplerate/2));

% filtered signals
emg_br_filt = filter(b, a, emg_br_rec);
emg_bi_filt = filter(b, a, emg_bi_rec);
emg_tlong_filt = filter(b, a, emg_tlong_rec);
emg_pd_filt = filter(b, a, emg_pd_rec);
emg_ad_filt = filter(b, a, emg_ad_rec);



figure
plot(emg_br_raw);
hold on
plot(emg_br_filt);
hold off
legend('BR');

figure
plot(emg_bi_raw);
hold on
plot(emg_bi_filt);
hold off
legend('BI');

figure
plot(emg_tlong_raw);
hold on
plot(emg_tlong_filt);
hold off
legend('TLONG');

figure
plot(emg_pd_raw);
hold on
plot(emg_pd_filt);
hold off
legend('PD');

figure
plot(emg_ad_raw);
hold on
plot(emg_ad_filt);
hold off
legend('AD');