
disp(gs_mps.AL);
N = 0;
delta = 1;
stag_h_field = 0;
h_field = 0;

charges = U1([1 -1]);
pspace = GradedSpace.new(charges, [1 1], false);
trivspace = GradedSpace.new(U1(0), 1, false);

%%

magn = get_magnetisation('XXZ', gs_mps, pspace, trivspace, 1.5, false, false);
fprintf('Magnetisation is %s \n', magn)

stag_magn = get_magnetisation('XXZ', gs_mps, pspace, trivspace, 1.5, false, true);
fprintf('Staggered magnetisation is %s \n', stag_magn)
% above converges to 1

stag_magn_new = get_magnetisation('XXZ', gs_mps, pspace, trivspace, 2, false, true);
fprintf('New staggered magnetisation is %s \n', stag_magn_new)
% above converges to 2

stag_magn = get_magnetisation('XXZ', gs_mps, pspace, trivspace, 2.1, false, true);
fprintf('Staggered magnetisation is %s \n', stag_magn)
