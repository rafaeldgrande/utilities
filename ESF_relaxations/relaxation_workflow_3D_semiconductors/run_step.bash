#!/bin/bash

export OMP_NUM_THREADS=1

module load qe/7.2
BGWBIN="/home1/08948/rdgrande/packages/BerkeleyGW-4.0/bin"
FORCES="/home1/08948/rdgrande/packages/excited_state_forces/excited_forces.py"
MYPYTHON="/home1/08948/rdgrande/miniconda3/bin/python"

# preparing QE and BGW files
bash ../scripts/00-create_directories.bash 
python ../scripts/01-generate_QE_inputs.py 
python ../scripts/02-generate_bgw_inputs.py
bash ../scripts/03-generate_links.bash

timestamp() { date +"%s"; }
startFULL=$(timestamp)

# --- 1-scf ---
cd 1-scf/
start=$(timestamp)
echo "[INFO] Starting 1-scf at $(date)"
ibrun -n 20 pw.x -pd .true. < qe.in &> qe.out
ibrun -n 1 pw2bgw.x < pw2bgw.in &> pw2bgw.out
$BGWBIN/wfn2hdf.x BIN wfn.complex WFN.h5
ibrun ph.x < ph.in &> ph.out
ibrun -n 2 dynmat.x < dynmat.in &> dynmat.log
end=$(timestamp)
echo "[INFO] Finished 1-scf at $(date) — Duration: $((end - start)) seconds"
cd ../

# --- 2-wfn_co ---
cd 2-wfn_co/
start=$(timestamp)
echo "[INFO] Starting 2-wfn_co at $(date)"
ibrun -n 20 pw.x -pd .true. < qe.in &> qe.out
ibrun -n 1 pw2bgw.x < pw2bgw.in &> pw2bgw.out
end=$(timestamp)
echo "[INFO] Finished 2-wfn_co at $(date) — Duration: $((end - start)) seconds"
cd ../

# --- 3-wfnq_co ---
cd 3-wfnq_co/
start=$(timestamp)
echo "[INFO] Starting 3-wfnq_co at $(date)"
ibrun -n 20 pw.x -pd .true. < qe.in &> qe.out
ibrun -n 1 pw2bgw.x < pw2bgw.in &> pw2bgw.out
$BGWBIN/wfn2hdf.x BIN wfn.complex WFN.h5
end=$(timestamp)
echo "[INFO] Finished 3-wfnq_co at $(date) — Duration: $((end - start)) seconds"
cd ../

# --- 4-parabands ---
cd 4-parabands/
start=$(timestamp)
echo "[INFO] Starting 4-parabands at $(date)"
ibrun $BGWBIN/parabands.cplx.x &> parabands.out
end=$(timestamp)
echo "[INFO] Finished 4-parabands at $(date) — Duration: $((end - start)) seconds"
cd ../

# --- 5-epsilon ---
cd 5-epsilon
start=$(timestamp)
echo "[INFO] Starting 5-epsilon at $(date)"
ibrun $BGWBIN/epsilon.cplx.x &> epsilon.out
end=$(timestamp)
echo "[INFO] Finished 5-epsilon at $(date) — Duration: $((end - start)) seconds"
cd ../

# --- 6-sigma ---
cd 6-sigma/
start=$(timestamp)
echo "[INFO] Starting 6-sigma at $(date)"
ibrun $BGWBIN/sigma.cplx.x &> sigma.out
end=$(timestamp)
echo "[INFO] Finished 6-sigma at $(date) — Duration: $((end - start)) seconds"
cd ../

# --- 7-kernel ---
cd 7-kernel/
start=$(timestamp)
echo "[INFO] Starting 7-kernel at $(date)"
ibrun $BGWBIN/kernel.cplx.x &> kernel.out
end=$(timestamp)
echo "[INFO] Finished 7-kernel at $(date) — Duration: $((end - start)) seconds"
cd ../

# --- 8-absorption ---
cd 8-absorption/
start=$(timestamp)
echo "[INFO] Starting 8-absorption at $(date)"
ibrun $BGWBIN/absorption.cplx.x &> absorption.out
end=$(timestamp)
echo "[INFO] Finished 8-absorption at $(date) — Duration: $((end - start)) seconds"
cd ../

# ---excited state forces and displacements calculation-----
start=$(timestamp)
echo "[INFO] Starting ESF at $(date)"
bash ../scripts/04-excited_forces.bash
end=$(timestamp)
echo "[INFO] Finished ESF at $(date) — Duration: $((end - start)) seconds"

# apply displacements
python ../scripts/apply_displacements.py

endFULL=$(timestamp)
echo "[INFO] Finished all calculations at $(date) — Duration: $((endFULL - startFULL)) seconds"

