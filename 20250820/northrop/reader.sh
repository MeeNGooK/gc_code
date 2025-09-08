#!/bin/bash

# ===== 1. NetCDF 모듈 로드 =====
module load netcdf-c/4.9.2
module load netcdf-fortran/4.6.1

# ===== 2. 빌드용 디렉토리 및 출력 파일 =====
MODDIR=build
OUTFILE=main.exe

mkdir -p "$MODDIR"

# ===== 3. NetCDF 관련 모듈 먼저 컴파일 =====
echo ">>> Compiling NetCDF reader..."
gfortran -c nc_reader.f90 -J "$MODDIR" -I$NETCDF_DIR/include

# ===== 4. 나머지 소스 파일 목록 (추가 소스도 쉽게 넣을 수 있음) =====
SRC=(
  modules/fourier1.f90
  main.f90

)

# ===== 5. 나머지 소스 컴파일 & 링크 (nc_reader.o 포함) =====
echo ">>> Compiling other sources..."
gfortran -I "$MODDIR" "${SRC[@]}" nc_reader.o -J "$MODDIR" \
    -I$NETCDF_DIR/include -L$NETCDF_DIR/lib -lnetcdff -o "$OUTFILE"

compile_status=$?

# ===== 6. 실행 =====
if [ $compile_status -eq 0 ]; then
  echo ">>> Compilation successful. Running $OUTFILE..."
  ./"$OUTFILE"
else
  echo ">>> Compilation failed. $OUTFILE will not run."
fi
