#!/bin/bash

export LC_ALL=C
sort --parallel 16 -S 20% -T . --batch-size=1021 $@
