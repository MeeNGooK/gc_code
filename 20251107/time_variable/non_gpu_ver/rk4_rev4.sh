#!/bin/bash

MODDIR=build
OUTFILE=main.exe

mkdir -p "$MODDIR"

# ===== 소스 파일 목록 =====
SRC=(
  set_values.cuf
  modules/trig_compute.cuf
  modules/spline2.cuf
  modules/fourier3.cuf
  modules/spline.cuf
  modules/fourier1.cuf
  reader/nc_reader.cuf
  rk4_revised4.cuf
)

echo ">>> Compiling with NVFORTRAN..."

nvfortran \
    -I"$MODDIR" -module "$MODDIR" \
    -I/home/rjrj524/modules/include \
    "${SRC[@]}" \
    -o "$OUTFILE" \
    $(/home/rjrj524/modules/bin/nf-config --flibs)

compile_status=$?

if [ $compile_status -eq 0 ]; then
  echo ">>> Compilation successful. Running $OUTFILE..."
  export LD_LIBRARY_PATH=/home/rjrj524/modules/lib:$LD_LIBRARY_PATH
  # base / step(interval)



  outfile="dat_e.dat"
  # u, v_perp, outfilename, timestep, iteration, verbose, m/m_e
  echo ">>> Running with arguments: 5.0e5 3.4e6 '$outfile' 1.0e-9 5000000 100 1"
  ./"$OUTFILE" 5.0e5 3.4e6 "$outfile" 1.0e-10 5000000 100 1

else
  echo ">>> Compilation failed. $OUTFILE will not run."

fi