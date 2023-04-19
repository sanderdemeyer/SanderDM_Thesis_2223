function mps = get_Hubbard_mps(P, Q, kwargs)
    arguments
        P
        Q
        kwargs.len = []
        kwargs.system = {'1D'}
        kwargs.D = 1
        kwargs.symmetries = 'U1_U1'
    end
    if strcmp(kwargs.symmetries, 'U1_U1')
        [pspace, vspaces, ~] = get_spaces_Hubbard_asymmetric(P, Q, 'D', kwargs.D);
    elseif strcmp(kwargs.symmetries, 'U1_SU2')
        [pspace, vspaces, ~] = get_spaces_Hubbard_SU2(P, Q, 'D', kwargs.D);
    elseif strcmp(kwargs.symmetries, 'None_SU2')
        [pspace, vspaces, ~] = get_spaces_Hubbard_None_SU2('D1', kwargs.D, 'D2', kwargs.D);
    else
        error('TBA');
    end
    if strcmp(kwargs.system{1}, 'Cylinder')
        N = kwargs.system{2};
        args = cell(2, max(N,2));
        for i = 1:max(N,2)
            args{1,i} = pspace;
            args{2,i} = vspaces((mod(i-1,length(vspaces))+1));
        end
        mps = UniformMps.randnc(args{:});
    elseif strcmp(kwargs.system{1}, 'Cylinder_multiple_rungs')
        N = kwargs.system{2};
        rungs = kwargs.system{3};
        args = cell(2, max(rungs*N,2));
        for i = 1:max(rungs*N,2)
            args{1,i} = pspace;
            args{2,i} = vspaces((mod(i-1,length(vspaces))+1));
        end
        mps = UniformMps.randnc(args{:});
    elseif strcmp(kwargs.system{1}, '1D') || strcmp(kwargs.system{1}, 'Helix')
        if isempty(kwargs.len)
            len = length(vspaces);
        else
            len = kwargs.len;
        end
        args = cell(2, length(vspaces));
        for i = 1:len
            args{1,i} = pspace;
            args{2,i} = vspaces((mod(i-1,length(vspaces))+1));
        end
        mps = UniformMps.randnc(args{:});
    else
        error('Invalid system: %s', kwargs.system{1});
    end
end