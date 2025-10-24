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
  ./"$OUTFILE"
else
  echo ">>> Compilation failed. $OUTFILE will not run."
fi
