function [d0,m0,normresult] = BalancedMultiLSAMF_new(d,m,N0,l,Norm)

% 均衡多道自适应匹配滤波
% 输入：
%   d m分别为实际数据以及预测的多次波（地表散射波），N0为每次计算时选取的道数(单侧)
%   l为自适应滤波器时窗长度, 
%   norm可以是1代表基于L1范数，2代表基于L2范数 
%   基于L1范数的算法参照 http://sepwww.stanford.edu/public/docs/sep103/antoine2/paper_html/node4.html#SECTION00022000000000000000
% 输出：
%   d0 m0分别为滤波后得到的有效波形以及多次波（地表散射波）
% d m d0 dm: nt * ntrace
% 相比于原来的修改：
%   为了使滤波器s能够使原始波形负向移动，将M改为了大小为n*(2l+1)的矩阵：
%   m(l+1)  m(l)    ...  m(1)  0     ...  0
%   m(l+2)  m(l+1)  ...  m(2)  m(1)  ...  0
%   .       .        .   .     .      .   .
%   .       .        .   .     .      .   m(1)
%   m(n)    .        .   .     .      .   .
%   0       m(n)     .   .     .      .   .
%   .       .        .   .     .      .   .
%   0       0       ...  m(n)  .     ...  m(n-l)
%   |--------l---------|  1  |------l------|
%             

Ncgls = 1e3;

[s1,s2] = size(d);
m0 = zeros(s1,s2);
if Norm == 2    %基于L2范数
    normresult = 12345;
    for itrace = 1:s2
        bgind = max(1,itrace - N0);
        edind = min(itrace + N0,s2);
        Mwhole = zeros(2*l+1);
        dwhole = zeros(2*l+1,1);
        for ik = bgind:edind
            Mi = zeros(s1,2*l+1);
            mik = m(:,ik);
            
            for is2 = 1:l+1
                Mi(1:s1-l+is2-1,is2) = mik(l+1-(is2-1):s1);
            end
            for is2 = l+2:2*l+1
                Mi(is2-l:s1,is2) = mik(1:l+s1-is2+1);
            end
            
            Mwhole = Mwhole + Mi' * Mi;
            dwhole = dwhole + Mi' * d(:,ik);

        end
        X = cgls(Mwhole,dwhole,Ncgls);
        if isnan(X(1,1))
            s = zeros(l,1);
        else
            s = X(:,end);
            while isnan(s(1))
                X = X(:,1:end-1);
                s = X(:,end);
            end
        end
        m0temp = conv(m(:,itrace),s);
        m0(:,itrace) = m0temp(1+l:s1+l);
        
    end
    d0 = d - m0;
    
elseif Norm == 1    %基于L1范数 需要迭代
    epsilon = max(abs(d(:)))/100;
    maxniter = 10;
    normresult = zeros(maxniter,s2);
    for itrace = 1:s2
        bgind = max(1,itrace - N0);
        edind = min(itrace + N0,s2);
        
        l1norm = sum(abs(d(:,bgind:edind)),[1 2]);
        l1normbf = l1norm + 1;
        
        niter = 1;
%         while l1norm < l1normbf && niter <= maxniter
        while niter <= maxniter
            Mwhole = zeros(2*l+1);
            dwhole = zeros(2*l+1,1);
            for ik = bgind:edind
                Mi = zeros(s1,2*l+1);
                mik = m(:,ik);

                for is2 = 1:l+1
                    Mi(1:s1-l+is2-1,is2) = mik(l+1-(is2-1):s1);
                end
                for is2 = l+2:2*l+1
                    Mi(is2-l:s1,is2) = mik(1:l+s1-is2+1);
                end
                if niter == 1
                    r = d(:,ik);
                else
                    r = d(:,ik) - Mi * s;
                end
                Wi = diag((1 + (r./epsilon).^2).^-0.25);
                Mwhole = Mwhole + Mi' * (Wi' * Wi) * Mi;
                dwhole = Mi' * (Wi' * Wi) * d(:,ik);
            end
            X = cgls(Mwhole,dwhole,Ncgls);
            if isnan(X(1,1))
                s = zeros(l,1);
            else
                s = X(:,end);
                while isnan(s(1))
                    X = X(:,1:end-1);
                    s = X(:,end);
                end
            end
            rM = zeros(s1,edind - bgind + 1);
            for ik = bgind:edind
                temp = conv(m(:,itrace),s);
                rM(:,ik-bgind+1) = d(:,ik) - temp(l+1:s1+l);
            end
            l1normbf = l1norm;
            l1norm = sum(abs(rM(:)));
            normresult(niter,itrace) = l1norm;
            niter = niter + 1;
        end
        
        m0temp = conv(m(:,itrace),s);
        m0(:,itrace) = m0temp(l+1:s1+l);
    end
    d0 = d - m0;
end