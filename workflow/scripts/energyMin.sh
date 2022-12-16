#!/bin/bash

set +e
obminimize -o pdbqt -cg -c 1e-5 -n 300 $1 | sed '/^TORSDOF/a ENDMDL' | sed '/^REMARK  Name =/i MODEL' > $2
if [$? -eq 1]
  then
    echo 'energyMin error on file: {input}'
    exit 0
else
  exit 0
fi
