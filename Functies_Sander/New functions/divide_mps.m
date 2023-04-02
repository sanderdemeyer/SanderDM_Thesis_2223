function mps_new = divide_mps(gs_mps, max)
    % makes an (N*rungs)-site mps based on a N-site mps
    args{4} = gs_mps.AC(1:max);
    args{3} = gs_mps.C(1:max);
    args{2} = gs_mps.AR(1:max);
    args{1} = gs_mps.AL(1:max);
    mps_new = UniformMps(args{:});
end