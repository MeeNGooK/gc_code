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
  base=20000
  step=2000

  for i in {1..10}; do
    arg2=$(echo "$base + ($i - 1) * $step" | bc -l)
    outfile="dat${i}.dat"
  # u, v_perp, outfilename, timestep, iteration, verbose, m/m_e
    echo ">>> Running with arguments: 5.0e5 $arg2 '$outfile' 1.0e-9 1000000 100 1"
    ./"$OUTFILE" $arg2 5.0e5 "$outfile" 1.0e-8 1000000 1 1836
  done

else
  echo ">>> Compilation failed. $OUTFILE will not run."
fi
