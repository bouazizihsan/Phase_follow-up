function x = constellation(be)
for i=1:length(be)
    switch be(i)
        case 0
            x(i) = -1+1j;
        case 1
            x(i) = -1-1j;
        case 2
            x(i) =  1 +1j;
        case 3
            x(i) = 1-1j;
    end
end
end