function [s_tensor] = tensorSignal(s,M)
N=length(s);
s=tensor(s);
s_t=s;
for ii=1:M-2
    s_t=squeeze(ttt(s_t,s));
end
otherdims = repmat({':'},1,ndims(s_t)-1);
s_tensor=zeros([N,1,repmat(N,1,ndims(s_t)-1)]);
s_tensor(:,1,otherdims{:})=double(s_t);
