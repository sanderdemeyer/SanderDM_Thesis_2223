function mps = get_Hubbard_mps_None_SU2(kwargs)
    arguments
        kwargs.system = {'1D'}
        kwargs.D1 = 1
        kwargs.D2 = 1
    end
    [pspace, vspaces, ~] = get_spaces_Hubbard_None_SU2('D1', kwargs.D1, 'D2', kwargs.D2);

    if strcmp(kwargs.system{1}, 'Cylinder')
        N = kwargs.system{2};
        args = cell(2, max(N,2));
        for i = 1:max(N,2)
            args{1,i} = pspace;
            args{2,i} = vspaces((mod(i-1,2)+1));
        end
        mps = UniformMps.randnc(args{:});
    elseif strcmp(kwargs.system{1}, 'DoubleCylinder')
        if P ~= 1 || Q ~= 1
            error('Only implemented for P = Q = 1');
        end
        N = kwargs.system{2};
        args = cell(2, max(2*N,2));
        for i = 1:max(2*N,2)
            args{1,i} = pspace;
            args{2,i} = vspaces((mod(i-1,2)+1));
        end
        mps = UniformMps.randnc(args{:});
    elseif strcmp(kwargs.system{1}, '1D')
        args = cell(2, length(vspaces));
        for i = 1:length(vspaces)
            args{1,i} = pspace;
            args{2,i} = vspaces(i);
        end
        mps = UniformMps.randnc(args{:});
    else
        error('System not implemented')
    end
end