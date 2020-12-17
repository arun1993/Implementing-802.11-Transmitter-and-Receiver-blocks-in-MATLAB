function [map_out]=demapper(data,mode,c)

m_out=data;
% mode : Modulation order in power of 2 (1/2/4/6 = BSPK/QPSK/16-QAM/64-QAM)
% c = 1;

switch mode
    case 1
        b=c*[1 -1];
    case 2
        b=c*[1+1i -1+1i 1-1i -1-1i];
    case 4
        b=c*[1+1i 1+3i 1-1i 1-3i 3+1i 3+3i 3-1i 3-3i -1+1i -1+3i -1-1i -1-3i -3+1i -3+3i -3-1i -3-3i];
    case 6
        b=c*[3+3i 3+1i 3+5i 3+7i 3-3i 3-1i 3-5i 3-7i 1+3i 1+1i 1+5i 1+7i 1-3i 1-1i 1-5i 1-7i 5+3i 5+1i 5+5i 5+7i 5-3i 5-1i 5-5i 5-7i 7+3i 7+1i 7+5i 7+7i 7-3i 7-1i 7-5i 7-7i -3+3i -3+1i -3+5i -3+7i -3-3i -3-1i -3-5i -3-7i -1+3i -1+1i -1+5i -1+7i -1-3i -1-1i -1-5i -1-7i -5+3i -5+1i -5+5i -5+7i -5-3i -5-1i -5-5i -5-7i -7+3i -7+1i -7+5i -7+7i -7-3i -7-1i -7-5i -7-7i];
    otherwise
        error('wrong choice');
end

temp = zeros(size(m_out));

for k=1:length(m_out)
    z=find( (abs(b-m_out(k))).^2-min(abs((b-m_out(k))).^2)==0);
    if length(z)==1
        temp(k)=z-1;
    else
        temp(k)=z(randi(1,1,[1,length(z)]))-1;
    end
end

dmod_out = de2bi(temp,mode);
map_out = reshape(dmod_out',1,mode*length(temp));

