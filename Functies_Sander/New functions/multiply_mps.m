function mps_new = multiply_mps(gs_mps, rungs)
    % makes an (N*rungs)-site mps based on a N-site mps
    args{4} = repmat(gs_mps.AC, 1, rungs);
    args{3} = repmat(gs_mps.C, 1, rungs);
    args{2} = repmat(gs_mps.AR, 1, rungs);
    args{1} = repmat(gs_mps.AL, 1, rungs);
    mps_new = UniformMps(args{:});
end