function m_hist = matrixify(v,n,m)
    % dim 1 of v should be the long dimension, like an ode45 ouput
    m_hist = reshape(v',n,m,[]);
end