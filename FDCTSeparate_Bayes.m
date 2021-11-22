function [d0,m0,d0total,m0total] = FDCTSeparate_Bayes(d,m,L,lambda1,lambda2,ita)

% Saab, R., D. Wang, O. Yilmaz, and F. J. Herrmann, 2007, Curvelet-based
% primary-multiple separation from a bayesian perspective: SEG
% International Exposition and 77th Annual Meeting, 2510?2514.
%
% d:  whole data
% m:  synthetic data to be separated from d
% d0: useful data remained
% m0: matching of synthetic data in real data
% size: 
%   d m d0 m0: nt * ntrace
% L:  max number of iterations
% lamda1 lambda2 ita defined as in the above paper

[M,N] = size(d);
S  = @(x) fdct_wrapping(x,1,1);
ST = @(x) ifdct_wrapping(x,1,M,N);

% initialization
x1 = S(d-m);
x2 = S(m);
b1 = d - m;     % estimation (no change during iteration)
b2 = m;         % estimation (no change during iteration)
m0 = m;
d0 = d - m;
nscale = size(x1,2);
w1 = S(zeros(size(m)));
w2 = w1;
for is = 1:nscale
    nangle = size(x1{is},2);
    for ia = 1:nangle
        w1{is}{ia} = lambda1 * abs(x2{is}{ia}) / (2*ita);
        w2{is}{ia} = lambda2 * abs(x1{is}{ia}) / (2*(1+ita));
    end
end

k = 0;

m0total = zeros(M,N,L);    % record the m0 of every outer loop
d0total = m0total;

% plot
% figure;
% subplot(2,2,1)
% imagesc(fliplr(d)./max(d(:)))
% caxis([-1 1])
% colormap redblue
% title('Real data')
% subplot(2,2,2)
% imagesc(fliplr(m)./max(d(:)))
% caxis([-1 1])
% colormap redblue
% title('1layer Synthetic data')

while k < L
    k = k + 1;
    disp(['Calculating the ',num2str(k),'th loop']);
    
    temp1 = b2 - m0 + b1 - d0;
    temp1 = CdomainSum(S(temp1),x1);
    x1 = Softthreshold(temp1,w1);
    d0 = ST(x1);
    
    temp2 = b2 - m0 + ita/(1+ita)*(b1 - d0);
    temp2 = CdomainSum(S(temp2),x2);
    x2 = Softthreshold(temp2,w2);
    m0 = ST(x2);
    
    m0total(:,:,k) = m0;
    d0total(:,:,k) = d0;
end

% subplot(2,2,3)
% imagesc(fliplr(d0)./max(d(:)))
% caxis([-1 1])
% colormap redblue
% title('Filtered Result')
% subplot(2,2,4)
% imagesc(fliplr(m0)./max(d(:)))
% caxis([-1 1])
% colormap redblue
% title('Matched Result')

end

% Calculate the sum of coefficients in curvelet domain
function s = CdomainSum(a,b)

s = a;
nscale = size(a,2);
for is = 1:nscale
    nangle = size(a{is},2);
    for ia = 1:nangle
        s{is}{ia} = s{is}{ia} + b{is}{ia};
    end
end

end

% Calculate Soft threshold
function x = Softthreshold(xk,w)

x = xk;
nscale = size(xk,2);
for is = 1:nscale
    nangle = size(xk{is},2);
    for ia = 1:nangle
        x{is}{ia} = sign(xk{is}{ia}) .* max(0,abs(xk{is}{ia}) - abs(w{is}{ia}));
    end
end

end