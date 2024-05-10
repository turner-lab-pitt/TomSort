function wavelet_coeffs = do_wavelet_decomposition(waveshapes)

if ~exist('wavedec','file') % Looks for Wavelets Toolbox
    error('extract_wave_features depends on the Wavelet toolbox for wavelet decomposition');
end
scales = 3;

num_spikes = size(waveshapes,1);

% [~,Hi_D]=wfilters('db12','d');
% V2=[-0.0000629061180000,0.0003436319050000,-0.0004539566200000,-0.0009448971360000,0.0028438345470000,0.0007081375040000,-0.0088391034090000,0.0031538470560000,0.0196872150100000,-0.0148534480050000,-0.0354703986070000,0.0387426192930000,0.0558925236910000,-0.0777097509020000,-0.0839288843660000,0.1319716614170000,0.1350842271290000,-0.1944504717660000,-0.2634948024880000,0.2016121617750000,0.6356010598720000,0.5727977932110000,0.2501841295050000,0.0457993341110000];
%wavelet features
% [c,~]=wavedec(feature_spikes(1,:),scales,'haar');
% [c,~]=wavedec(feature_spikes(1,:),scales,V2,Hi_D);
[c,~]=wavedec(waveshapes(1,:,1),scales,'db6');
ls=length(c);
wavelet_coeffs=zeros(num_spikes,ls);
parfor i=1:num_spikes                                % Wavelet decomposition
%     [c,~]=wavedec(feature_spikes(i,:),scales,'haar');
%     [c,~]=wavedec(feature_spikes(i,:),scales,V2,Hi_D);
    [c,~]=wavedec(waveshapes(i,:),scales,'db6');
%     cc(i,1:ls)=c(1:ls);
    wavelet_coeffs(i,:)=c;
end
% output.features.wavelet_coeffs=wavelet_coeffs;

end