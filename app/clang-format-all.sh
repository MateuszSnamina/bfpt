#!/bin/bash


#find . \( -name '*.hpp' -or -name '*.cpp' \) -print0

find . \( -name '*.hpp' -or -name '*.cpp' \) -print0 | \
xargs -0 clang-format -i
