#!/bin/bash

set -xe

rm -f main
gcc -o main main.c -lraylib -lm -lpthread -Wall -Wextra
./main

