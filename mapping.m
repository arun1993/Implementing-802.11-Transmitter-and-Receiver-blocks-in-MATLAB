function [map_out]=mapping(data,mode,scale)

% mode : Modulation order in power of 2 (1/2/4/6 = BSPK/QPSK/16-QAM/64-QAM)
% scale : scaling up or down the modulation, default value = 1

input_seq = data;

switch mode
    case 1
        b=scale*[1 -1];
    case 2
        b=scale*[1+1i -1+1i 1-1i -1-1i];
    case 4
        b=scale*[1+1i 1+3i 1-1i 1-3i 3+1i 3+3i 3-1i 3-3i -1+1i -1+3i -1-1i -1-3i -3+1i -3+3i -3-1i -3-3i];
    case 6
        b=scale*[3+3i 3+1i 3+5i 3+7i 3-3i 3-1i 3-5i 3-7i 1+3i 1+1i 1+5i 1+7i 1-3i 1-1i 1-5i 1-7i 5+3i 5+1i 5+5i 5+7i 5-3i 5-1i 5-5i 5-7i 7+3i 7+1i 7+5i 7+7i 7-3i 7-1i 7-5i 7-7i -3+3i -3+1i -3+5i -3+7i -3-3i -3-1i -3-5i -3-7i -1+3i -1+1i -1+5i -1+7i -1-3i -1-1i -1-5i -1-7i -5+3i -5+1i -5+5i -5+7i -5-3i -5-1i -5-5i -5-7i -7+3i -7+1i -7+5i -7+7i -7-3i -7-1i -7-5i -7-7i];
    otherwise
        error('wrong choice');
end

count=1;
map_out = zeros(1,ceil(length(input_seq)/mode));

for i=1:(ceil(length(input_seq)/mode))
    temp=0;
    for j=1:mode
        temp=bitor(temp,bitshift(input_seq(count),(j-1)));
        count=count+1;
        if(count>length(input_seq))
            break;
        end
    end
    map_out(i)=b(temp+1);
end

