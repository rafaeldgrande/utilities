
prefix='si'

# 2-wfn_co
ln -sf ../../1-scf/$prefix.save/data-file-schema.xml 2-wfn_co/$prefix.save/
ln -sf ../../1-scf/$prefix.save/charge-density.dat 2-wfn_co/$prefix.save/
ln -sf ../../1-scf/$prefix.save/charge-density.hdf5 2-wfn_co/$prefix.save/

# 3-wfnq_co
ln -sf ../../1-scf/$prefix.save/data-file-schema.xml 3-wfnq_co/$prefix.save/
ln -sf ../../1-scf/$prefix.save/charge-density.dat 3-wfnq_co/$prefix.save/
ln -sf ../../1-scf/$prefix.save/charge-density.hdf5 3-wfnq_co/$prefix.save/

# 4-parabands
ln -sf ../2-wfn_co/wfn.complex 4-parabands/
ln -sf ../2-wfn_co/VKB 4-parabands/
ln -sf ../2-wfn_co/VSC 4-parabands/

# 5-epsilon
ln -sf ../4-parabands/WFN.h5 5-epsilon/WFN.h5
ln -sf ../3-wfnq_co/WFN.h5 5-epsilon/WFNq.h5

# 6-sigma
ln -sf ../2-wfn_co/RHO 6-sigma/RHO
ln -sf ../2-wfn_co/VXC 6-sigma/VXC
ln -sf ../4-parabands/WFN.h5 6-sigma/WFN_inner.h5 
ln -sf ../4-parabands/WFN.h5 6-sigma/WFN_outer.h5
ln -sf ../5-epsilon/eps0mat.h5 6-sigma/eps0mat.h5
ln -sf ../5-epsilon/epsmat.h5 6-sigma/epsmat.h5

# 7-kernel
ln -sf ../4-parabands/WFN.h5 7-kernel/WFN_co.h5
ln -sf ../5-epsilon/eps0mat.h5 7-kernel/eps0mat.h5
ln -sf ../5-epsilon/epsmat.h5 7-kernel/epsmat.h5

# 8-absorption
ln -sf ../4-parabands/WFN.h5 8-absorption/WFN_co.h5
ln -sf ../1-scf/WFN.h5 8-absorption/WFN_fi.h5
ln -sf ../5-epsilon/eps0mat.h5 8-absorption/eps0mat.h5
ln -sf ../5-epsilon/epsmat.h5 8-absorption/epsmat.h5
ln -sf ../6-sigma/eqp1.dat 8-absorption/eqp_co.dat
ln -sf ../7-kernel/bsemat.h5 8-absorption/bsemat.h5

# 9-excited-state-forces 
ln -sf ../8-absorption/eigenvectors.h5 9-excited-state-forces/
ln -sf ../8-absorption/eqp.dat 9-excited-state-forces/
ln -sf ../1-scf/_ph0/$prefix.phsave 9-excited-state-forces/elph_dir


