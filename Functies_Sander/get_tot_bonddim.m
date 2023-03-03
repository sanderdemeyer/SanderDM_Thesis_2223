function som = get_tot_bonddim(mps)
    som = 0;
    for i = 1:period(mps)
        som = som + mps.AL(i).var.dims(1);
    end
