% Opties: 
% 0 e- aanwezig: U1(0) en SU2(1)
% 1 e- aanwezig: U1(1) en SU2(2)
% 2 e- aanwezig: U1(2) en SU2(1)
% In totaal zijn er zo 4 opties:
% 0, up, down, en dubbel bezet.

a = ProductCharge(U1(0), SU2(1), fZ2(0));
b = ProductCharge(U1(1), SU2(2), fZ2(1));
c = ProductCharge(U1(2), SU2(1), fZ2(0));
charges = [a b c];

fusion_trees = FusionTree.new([2 2], charges, false, charges, false, charges, false, charges, false);
pspace = GradedSpace.new(charges, [1 1 1], false);
trivspace = GradedSpace.new(a, 1, false);

%%

% For atomic separation = 2å of 1D hydrogen chain: 
% t = 1.5 en abs(U/t) = 6
t = 1.5;
U = 9;
mu = 1000;

H = get_hamiltonian('Hubbard', fusion_trees, pspace, t, U, mu);
mpo = get_mpo(H, 0);

H1 = get_JMpo(mpo, 1);
%H1 = get_JMpo({mpo mpo},1) %uncomment for two-site operator

D = 100;
vspace = GradedSpace.new(charges, [D D D], false);

mps = UniformMps.randnc(pspace', vspace);
%mps = UniformMps.randnc([pspace' pspace'], [vspace vspace]); %uncomment
%for two-site operator



%%

alg = Vumps('which', 'smallestreal', 'miniter', 1, 'maxiter', 25, 'verbosity', Verbosity.iter);
[gs_mps, gs_energy] = fixedpoint(alg, H1, mps);

%%

numbers = get_numbers(gs_mps, pspace, trivspace, 1);

