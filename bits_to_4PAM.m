function symbols = bits_to_4PAM(bits)
    % Convert a stream of binary bits to 4PAM symbols
    % 00 -> +3
    % 01 -> +1
    % 11 -> -1
    
    
    % Initialize symbols vector
    symbols = zeros(length(bits)/2, 1);
    
    % Convert binary bits to 4PAM symbols
    for i = 1:length(bits)/2
        if bits(2*i-1) == 0 && bits(2*i) == 0
            symbols(i) = 3;
        elseif bits(2*i-1) == 0 && bits(2*i) == 1
            symbols(i) = 1;
        elseif bits(2*i-1) == 1 && bits(2*i) == 1
            symbols(i) = -1;

        end
    end
end