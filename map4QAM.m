function vectS = map4QAM(vectB)
L = length(vectB)
re =reshape(vectB,[2,L/2])';
be =bi2de(re, 'left-msb');

for i=1:length(be)
    vectS(i) = mapping(be(i));
end

end

function  z = mapping(be)
switch be
    case 0
        z = -1+1j;
    case 1
        z=  -1-1j;
    case 2
        z = 1 +1j;
    case 3
        z = 1-1j;
end

end