#!/bin/bash

# 빌드용 디렉토리 (모듈 .mod 파일 저장)
module load netcdf-c/4.9.2
module load netcdf-fortran/4.6.1

MODDIR=build

# 출력 파일 이름
OUTFILE=main.exe

# 소스 파일 목록
SRC=(
  modules/cubic_eval_kernels.cuf
  modules/spline2.cuf

  modules/fourier2.cuf
  modules/spline.cuf
  modules/fourier1.cuf
  reader/nc_reader.cuf
  rk4_revised.cuf
)

# build 디렉토리 없으면 생성
mkdir -p "$MODDIR"

# 컴파일
echo ">>> Compiling..."
nvfortran -module "$MODDIR" -I "$MODDIR" "${SRC[@]}" -o "$OUTFILE"
compile_status=$?

# 컴파일 성공 여부 확인
if [ $compile_status -eq 0 ]; then
  echo ">>> Compilation successful. Running $OUTFILE..."
  ./"$OUTFILE"
else
  echo ">>> Compilation failed. $OUTFILE will not run."
fi
