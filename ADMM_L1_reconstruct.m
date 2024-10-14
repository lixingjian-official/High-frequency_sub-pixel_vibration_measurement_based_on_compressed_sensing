function x=ADMM_L1_reconstruct(A,b)
    delta = 100;
    lambda = 10;
    iteratMax =100;
    [~,N]=size(A);
    e = ones(N,1);
    D_ = spdiags([e -e], 0:1, N,N);
    I = eye(N);
    d=D_*ones(N,1);
    p=ones(N,1)/delta;
    invDD=(A'*A+delta*I);
    for ii=1:iteratMax
        x=invDD\(A'*b+delta*(d-p));
        d=wthresh(x+p,'s',lambda/delta);
        p=p+x-d;
    end
end
