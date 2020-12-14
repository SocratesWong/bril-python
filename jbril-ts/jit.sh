#!/bin/bash
PATH=$PATH:`yarn global bin` bril2json < $1 | PATH=$PATH:`yarn global bin`  jbrili $2 > temp.txt
echo '{"op": "speculate"},'>temp2.txt
sed 's/$/,/' temp.txt>> temp2.txt 
echo '{"op": "commit"},'>>temp2.txt
echo '{"op": "ret"},'>>temp2.txt
echo '{"label": "fail"},'>>temp2.txt
PATH=$PATH:`yarn global bin` bril2json < $1 > json.txt
sed '/"instrs": \[/ r temp2.txt' json.txt > out.txt

