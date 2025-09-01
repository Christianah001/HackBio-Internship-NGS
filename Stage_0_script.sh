#!/bin/bash

 #1.Print your name
echo Funmilayo


#2.Create a folder titled your name
mkdir Funmilayo


#3. Create another new directory titled biocomputing and change to that directory with one line of command
 mkdir biocomputing; cd biocomputing

#4.Download these 3 files:
https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk
https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk

wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk -O wildtype.gbk
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk -O wildtype.gbk.1



#5. Move the .fna file to the folder titled your name
mv *.fna ../Funmilayo


#6. Delete the duplicate gbk file
rm *.gbk.1


 #7. Confirm if the .fna file is mutant or wild type (tatatata vs tata)
#8. If mutant, print all matching lines into a new file

#Define the path to the .fna file
fna_file="../Funmilayo/wildtype.fna"
# Define the output file path in the Funmilayo directory
output_file="../Funmilayo/mutant_sequence.txt"

# Check if the .fna file contains "tatatata"
if grep -q "tatatata" "$fna_file"; then
    echo "Mutant"
    # Print all matching lines containing "tatatata" into a new file
    grep "tatatata" "$fna_file" > "$output_file"
elif grep -q "tata" "$fna_file"; then
    echo "Wild Type"
else
    echo "No matching sequence found."
fi

ls ../Funmilayo/

head -n 5 ../Funmilayo/mutant_sequence.txt 


#9. Count number of lines (excluding header) in the .gbk file
wc -l *.gbk
sed '1d' *.gbk | wc -l

#10. Print the sequence length of the .gbk file. (Use the LOCUS tag in the first line)
awk '/LOCUS/ {print $3}' *.gbk



#11. Print the source organism of the .gbk file. (Use the SOURCE tag in the first line)
awk '/SOURCE/ {print $2, $3}' *.gbk


#12. List all the gene names of the .gbk file. Hint {grep '/gene='}
grep '/gene=' *.gbk
                    
