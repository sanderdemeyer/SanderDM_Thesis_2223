function [list_D, list_angles, list_charges, list_charges_strings] = plot_transfereigs_None_SU2(gs_mps, kwargs)
    arguments
        gs_mps
        kwargs.sectors = 20
        kwargs.howmany = 25
        kwargs.symmetries = 'None_SU2'
        kwargs.charges_old = {[1 2 1] [-23 1 0] [-23 3 0] [-15 2 1] [-15 4 1] [-7 1 0] [-7 3 0] [-7 5 0] [1 4 1] [1 6 1] [9 1 0] [9 3 0] [9 5 0] [17 2 1] [17 4 1] [17 6 1] [25 1 0] [25 3 0] [25 5 0] [33 2 1]}
        kwargs.charges = {[-3 2 1] [-3 4 1] [-2 1 0] [-2 3 0] [-2 5 0] [-1 2 1] [-1 4 1] [-1 6 1] [-1 8 1] [0 1 0] [0 3 0] [0 5 0] [0 7 0] [1 2 1] [1 4 1]}
    end

    if length(kwargs.howmany) == 1
        kwargs.howmany = repmat(kwargs.howmany, 1, kwargs.sectors);
    else
        assert(length(kwargs.howmany) == kwargs.sectors, ...
            'Length of howmany should be the same as the number of sectors');
    end

    list_D = cell(1, kwargs.sectors);
    list_angles = cell(1, kwargs.sectors);
    list_angles_norm = cell(1, kwargs.sectors);
    list_charges = cell(1, kwargs.sectors);
    list_charges_strings = cell(1, kwargs.sectors);
    figure;
    for i = 1:kwargs.sectors
        disp(i);
        if strcmp(kwargs.symmetries, 'None_SU2')
            charge = ProductCharge(SU2(i), fZ2(mod(i+1,2)));
        else
            charges = kwargs.charges{i};
            charge = ProductCharge(U1(charges(1)), SU2(charges(2)), fZ2(charges(3)));
        end
        list_charges{i} = charge;
        [~, D] = transfereigs(gs_mps, gs_mps, kwargs.howmany(i), 'largestabs', 'Charge', charge); %, 'Type', 'l_LL');
        Diag = diag(D);
        list_D{i} = Diag;
        list_angles{i} = angle(Diag);
        list_angles_norm{i} = angle(Diag)/pi;
        list_charges_strings{i} = 'SU2(' + string(i) + ') x fZ2(' + string(mod(i+1,2)) + ')';
        scatter(real(Diag), imag(Diag));
        hold on
    end
    legend(list_charges_strings{:});
    viscircles([0 0], 1, 'Color', 'black');
    hold on;
    x = linspace(-1, 1, 500);
    for i = 0:7
        plot(x*cos(i*pi/8), x*sin(i*pi/8), 'Color', 'black');
        hold on;
    end
    axis equal;
    hold off;
    
    figure;
    histogram(cat(1, list_angles_norm{:}), 300);
end