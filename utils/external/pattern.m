function output = pattern(B, A)
    SIZE = length(B) - length(A);
    match = zeros(1, SIZE);
    for i=1:SIZE
        match(i) = all(B(i:i-1+length(A)) == A);
    end
    output = find(match == 1);
end