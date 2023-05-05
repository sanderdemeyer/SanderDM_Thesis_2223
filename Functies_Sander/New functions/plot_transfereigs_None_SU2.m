function [list_D, list_angles, list_charges, list_charges_strings] = plot_transfereigs_None_SU2(gs_mps, kwargs)
    arguments
        gs_mps
        kwargs.sectors = 5
        kwargs.howmany = 25
    end

    if length(kwargs.howmany) == 1
        kwargs.howmany = repmat(kwargs.howmany, 1, kwargs.sectors);
    else
        assert(length(kwargs.howmany) == kwargs.sectors, ...
            'Length of howmany should be the same as the number of sectors');
    end

    list_D = cell(1, kwargs.sectors);
    list_angles = cell(1, kwargs.sectors);
    list_charges = cell(1, kwargs.sectors);
    list_charges_strings = cell(1, kwargs.sectors);
    figure;
    for i = 1:kwargs.sectors
        disp(i);
        charge = ProductCharge(SU2(i), fZ2(mod(i+1,2)));
        list_charges{i} = charge;
        [~, D] = transfereigs(gs_mps, gs_mps, kwargs.howmany(i), 'largestabs', 'Charge', charge);
        Diag = diag(D);
        list_D{i} = Diag;
        list_angles{i} = angle(Diag);
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
    histogram(cat(1, list_angles{:}), 300);
end