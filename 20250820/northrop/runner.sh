#!/bin/bash

# 빌드용 디렉토리 (모듈 .mod 파일 저장)
MODDIR=build

# 출력 파일 이름
OUTFILE=main.exe

# 소스 파일 목록
SRC=(
  modules/crossp.cuf
  modules/electric_field.cuf
  modules/magnetic_field.cuf
  modules/deriv.cuf
  modules/constants.cuf
  modules/rk4.cuf
  single_lorentz.cuf
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
