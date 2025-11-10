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
  echo ">>> Running with arguments: 5.0e6 5.0e5 'dat1.dat'" 1.0e-8 1000000 100
  ./"$OUTFILE" 1.0e6 2.0e7 'dat1.dat' 1.0e-9 5000000 100
  echo ">>> Running with arguments: 5.0e6 6.0e5 'dat2.dat'" 1.0e-9 1000000 100
  ./"$OUTFILE" 1.0e6 2.5e7 'dat2.dat' 1.0e-9 5000000 100
  echo ">>> Running with arguments: 5.0e6 5.0e5 'dat3.dat'" 1.0e-9 1000000 100
  ./"$OUTFILE" 1.0e6 3.0e7 'dat3.dat' 1.0e-9 5000000 100
  echo ">>> Running with arguments: 5.0e6 6.0e55 'dat4.dat'" 1.0e-9 1000000 100
  ./"$OUTFILE" 1.0e6 3.5e7 'dat4.dat' 1.0e-9 5000000 100
  echo ">>> Running with arguments: 5.0e6 5.0e5 'dat5.dat'" 1.0e-9 1000000 100
  ./"$OUTFILE" 1.0e6 4.0e7 'dat5.dat' 1.0e-9 5000000 100
  echo ">>> Running with arguments: 5.0e6 6.0e5 'dat6.dat'" 1.0e-9 1000000 100
  ./"$OUTFILE" 1.0e6 4.5e7 'dat6.dat' 1.0e-9 5000000 100

else
  echo ">>> Compilation failed. $OUTFILE will not run."
fi
