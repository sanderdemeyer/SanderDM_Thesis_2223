function mps = get_Hubbard_mps(P, Q, kwargs)
    arguments
        P
        Q
        kwargs.system = {'1D'}
        kwargs.D = 1
    end
    %if P == 153488856512544475
    %    [pspace, vspaces, ~] = get_spaces_Hubbard_symmetric(P, Q, 'D1', kwargs.D, 'D2', kwargs.D);
    %else
    [pspace, vspaces, ~] = get_spaces_Hubbard_asymmetric(P, Q, 'D', kwargs.D);

    if strcmp(kwargs.system{1}, 'Cylinder')
        if P ~= 1 || Q ~= 1
            error('Only implemented for P = Q = 1');
        end
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
    else
        args = cell(2, length(vspaces));
        for i = 1:length(vspaces)
            args{1,i} = pspace;
            args{2,i} = vspaces(i);
        end
        mps = UniformMps.randnc(args{:});
    end
end
