#!/bin/bash -l

for FileName in *.pep
do
nentferner.pl -in=$FileName -out=$FileName.nent
sed -i 's/ .\+//g' $FileName.nent
sed -i 's/_c_s//g' $FileName.nent
done

rename .pep.nent .fa *.pep.nent

mkdir original_pep_files

mv *.pep ./original_pep_files/