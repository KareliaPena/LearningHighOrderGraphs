function [s_tensor] = polySignal(s,M)
N=length(s);
if M==3
    s_tensor=zeros(N,1,N);
    for i =1:N
        s_tensor(:,1,i)=s.^(i);
    end
end

