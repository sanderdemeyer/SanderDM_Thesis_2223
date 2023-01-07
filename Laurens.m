function [u1,n1,n2,hop]=FermionicOperators(choice)

    jk=1; dk=1; ak=+1;
    
    a_up=TensorSymm.Zeros({'Z2'},{dk [2;2] [2;2]},{jk [0;1] [0;1]},[-ak -1 +1],1);
    a_up.var{ismember(a_up.tree,[jk 0 1],'rows')}=[1 0; 0 0];
    a_up.var{ismember(a_up.tree,[jk 1 0],'rows')}=[0 0; 0 -1];
    c_up=TensorSymm.Zeros({'Z2'},{[2;2] [2;2] dk},{[0;1] [0;1] jk},[-1 +1 +ak],1);
    c_up.var{ismember(c_up.tree,[0 1 jk],'rows')}=[0 0; 0 -1];
    c_up.var{ismember(c_up.tree,[1 0 jk],'rows')}=[1 0; 0 0];
    a_do=TensorSymm.Zeros({'Z2'},{dk [2;2] [2;2]},{jk [0;1] [0;1]},[-ak -1 +1],1);
    a_do.var{ismember(a_do.tree,[jk 0 1],'rows')}=[0 1; 0 0];
    a_do.var{ismember(a_do.tree,[jk 1 0],'rows')}=[0 1; 0 0];
    c_do=TensorSymm.Zeros({'Z2'},{[2;2] [2;2] dk},{[0;1] [0;1] jk},[-1 +1 +ak],1);
    c_do.var{ismember(c_do.tree,[0 1 1],'rows')}=[0 0; 1 0];
    c_do.var{ismember(c_do.tree,[1 0 1],'rows')}=[0 0; 1 0];
    
    unit=TensorSymm.Zeros({'Z2'},{[2;2] [2;2]},{[0;1] [0;1]},[-1 +1],1);
    unit.var{1}=eye(2); unit.var{2}=eye(2);
    
    c_up=Fermionize(c_up); a_up=Fermionize(a_up);
    c_do=Fermionize(c_do); a_do=Fermionize(a_do);
    unit=Fermionize(unit);
    n_do=Contract({c_do,a_do},{[-1,1,2],[2,1,-2]});
    n_up=Contract({c_up,a_up},{[-1,1,2],[2,1,-2]});
    n_total=n_up+n_do;
    hop_do=Contract({c_do,a_do},{[-1,-3,1],[1,-2,-4]});
    hop_up=Contract({c_up,a_up},{[-1,-3,1],[1,-2,-4]});
    hop_total=hop_do+hop_up;
    
    switch choice
        case 'Z2'
        Dj_phys=[0; 1]; b_phys=[2;2]; sym={'Z2'};
        case 'su(2)xu(1)'
        Dj_phys=[0 1 0; 0 1 2; 1 2 1]; b_phys=[1;1;1]; sym={'Z2' 'su(2)' 'u(1)'};
        case 'su(2)'
        Dj_phys=[0 1; 1 2]; b_phys=[2;1]; sym={'Z2' 'su(2)'};
    end
    
    E1=TensorSymm.Zeros(sym,{b_phys b_phys},{Dj_phys Dj_phys},[-1 +1],1);
    E2=TensorSymm.Zeros(sym,{b_phys b_phys b_phys b_phys},{Dj_phys Dj_phys Dj_phys Dj_phys},[-1 -1 +1 +1],1);
    u1=Fermionize(FromArray(E1,double(unit)));
    n1=Fermionize(FromArray(E1,double(n_total)));
    hop=Fermionize(FromArray(E2,double(hop_total)));
    n2=Contract({n_total,unit},{[-1,-3],[-2,-4]})+Contract({unit,n_total},{[-1,-3],[-2,-4]});
    
end
    
    function A=Fermionize(A)
    A.group(1)=SymmetryGroup.MakeData('fermion');
    A.fermion=1;
end

