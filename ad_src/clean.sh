#!/bin/sh
# Removes the out-of-source build/ directory plus any stray in-source
# CMake files (leftover from running `cmake .` directly in ad_src/).
set -e
cd "$(dirname "$0")"

rm -rf build
rm -rf CMakeCache.txt CMakeFiles Makefile cmake_install.cmake

echo "ad_src cleaned."
