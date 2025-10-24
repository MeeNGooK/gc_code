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
  rk4_revised3.cuf
)

echo ">>> Compiling with NVFORTRAN..."

nvfortran \
    -I"$MODDIR" -module "$MODDIR" \
    -I/home/rjrj524/modules/include \
    "${SRC[@]}" \
    -o "$OUTFILE" \
    $(/home/rjrj524/modules/bin/nf-config --flibs)

compile_status=$?

# ===== 실행 =====
if [ $compile_status -eq 0 ]; then
  echo ">>> Compilation successful. Running $OUTFILE..."
  export LD_LIBRARY_PATH=/home/rjrj524/modules/lib:$LD_LIBRARY_PATH
  echo ">>> Running with arguments: 5.0e6 1.0e6 'dat11.dat'"
  ./"$OUTFILE" 5.0e5 1.40e6 'dat31.dat'
  echo ">>> Running with arguments: 5.0e6 2.0e6 'dat12.dat'"
  ./"$OUTFILE" 5.0e5 1.45e6 'dat32.dat'
  echo ">>> Running with arguments: 5.0e6 3.0e6 'dat13.dat'"
  ./"$OUTFILE" 5.0e5 1.50e6 'dat33.dat'
  echo ">>> Running with arguments: 5.0e6 4.0e5 'dat14.dat'"
  ./"$OUTFILE" 5.0e5 1.55e6 'dat34.dat'
  echo ">>> Running with arguments: 5.0e6 5.0e5 'dat15.dat'"
  ./"$OUTFILE" 5.0e5 1.60e6 'dat35.dat'
  echo ">>> Running with arguments: 5.0e6 6.0e5 'dat16.dat'"
  ./"$OUTFILE" 5.0e5 1.65e6 'dat36.dat'
  echo ">>> Running with arguments: 5.0e6 7.0e5 'dat17.dat'"
  ./"$OUTFILE" 5.0e5 1.70e6 'dat37.dat'
  echo ">>> Running with arguments: 5.0e6 8.0e6 'dat18.dat'"
  ./"$OUTFILE" 5.0e5 1.75e6 'dat38.dat'
  echo ">>> Running with arguments: 5.0e6 9.0e6 'dat39.dat'"
  ./"$OUTFILE" 5.0e5 1.80e6 'dat39.dat'
  echo ">>> Running with arguments: 5.0e6 1.0e7 'dat40.dat'"
  ./"$OUTFILE" 5.0e5 1.85e6 'dat40.dat'
else
  echo ">>> Compilation failed. $OUTFILE will not run."
fi
