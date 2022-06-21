function bsi = singleBSI(left, right)
    arguments
        left  (:,1) double
        right (:,1) double
    end
    % Check inputs:
    assert(length(left)==length(right));
    M = length(left);

    % Define Functions:
    Power = @(sig,i,j) abs(sig(j));
    R = @(i,j) Power(right,i,j);
    L = @(i,j) Power(left ,i,j);

    % Iterate
    sum_ = 0.00;
    for j = 1 : M
        Rij = R(0,j);
        Lij = L(0,j);
        val_ = ( Rij-Lij ) / ( Rij+Lij );
        sum_ = sum_ + abs(val_);
    end
    bsi = sum_ / M;    

end