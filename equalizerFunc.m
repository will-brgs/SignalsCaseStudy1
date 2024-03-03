function audio_output = equalizerFunc (input, gains, center_band, k_cut)
    Fs = 44.1e3;
    t = 0:1/Fs:0.25; %Sample timepoint vector

    R_Hi = zeros(5, 1);
    R_Lo = zeros(5, 1);
    for i = 1:5 % Calculate and allocate A and B values for Hi and Lo Coefficients
        % Calc Hipass R
        cutoff_Hi = center_band(i) - (k_cut * center_band(i));
        R_Hi(i) = 1/(2 * pi * C * cutoff_Hi);
        
        % Calc Lopass R
        cutoff_Lo = center_band(i) + (k_cut * center_band(i));
        R_Lo(i) = 1/(2 * pi * C * cutoff_Lo);
    
    end
    
    %Lowpass Filter Coefficients
    a_Lo = zeros(5, 2);
    a_Lo(:,1) = 1;
    a_Lo(:,2) = 1./(C.*R_Lo);
    
    b_Lo = 1./(C.*R_Lo);
    
    %HighPass Filter Coefficients
    a_Hi = zeros(5, 1);
    a_Hi(:,1) = 1;
    a_Hi(:,2) = 1./(C.*R_Hi);
    
    b_Hi = zeros(5, 2);
    b_Hi(:,1) = 1;

    audio_output = zeros(length(t), 1);
    
    for j = 1:5
        input_lo = lsim(b_Lo(j,:),a_Lo(j, :), input, t);
        input_band = lsim(b_Hi(j,:),a_Hi(j,:), input_lo, t);
        audio_output = audio_output + gains(j)*input_band;
    end      
end