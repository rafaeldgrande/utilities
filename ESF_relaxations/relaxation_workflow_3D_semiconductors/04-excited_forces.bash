
# necessary to not mess up my conda installation with frontera paths
unset PYTHONPATH PYTHONHOME LD_LIBRARY_PATH

echo -e "exc 1
jexc 1
eqp_file         eqp.dat
exciton_file     eigenvectors.h5
el_ph_dir        elph_dir/" > 9-excited-state-forces/forces.inp

MYPYTHON="/home1/08948/rdgrande/miniconda3/bin/python"
FORCESPATH="/home1/08948/rdgrande/packages/excited_state_forces"

cd 9-excited-state-forces
$MYPYTHON $FORCESPATH/excited_forces.py &> excited_state_forces.out 
cd ../

cd 10-displacements/

echo -e "file_out_QE ../1-scf/qe.out 
excited_state_forces_file ../9-excited-state-forces/forces_cart.out_1_1
eigvecs_file ../1-scf/eigvecs
qe_input ../1-scf/qe.in
limit_disp_eigvec_basis 0.3" > harmonic_approx.inp
$MYPYTHON $FORCESPATH/harmonic_extrapolation.py 
cd ../

