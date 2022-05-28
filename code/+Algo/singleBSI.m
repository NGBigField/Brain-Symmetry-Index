function bsi = singleBSI(L, R)
    arguments
        L (:,1) double
        R (:,1) double
    end
    % Check inputs:
    assert(length(L)==length(R));
    M = length(L);

    sum_ = 0.00;
    for j = 1 : M
        val_ = ( R(j)-L(j) ) / ( R(j)+L(j) );
        sum_ = sum_ + abs(val_);
    end
    bsi = sum_ / M;    

end