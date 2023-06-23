#!/bin/bash
echo Find Selection in all mammals
echo Required a fasta
echo To run use this script as ./Selection.sh fasta file
echo Checking the inputs 
cp $1 Fasta.fasta
#cut -d '_' -f1 Fasta.fasta >Fasta_H.fasta
cut -d '|' -f1 Fasta.fasta >Fasta_H.fasta
sed -i 's/_lcl//g' Fasta_H.fasta
rm Fasta.fasta
#echo Alignmnet
#java -jar macse_v2.04.jar -prog alignSequences -seq Fasta_H.fasta
#rm Fasta_H_AA.fasta
#mv Fasta_H_NT.fasta Fasta_AL.fasta
#rm Fasta_H.fasta
echo Remove internal stop condon
java -jar macse_v2.04.jar -prog exportAlignment -align Fasta_H.fasta -codonForFinalStop --- -codonForInternalStop NNN
rm Fasta_H_AA.fasta
mv Fasta_H_NT.fasta Fasta_AL_SR.fasta
rm Fasta_H.fasta
sed 's/\!/-/g' Fasta_AL_SR.fasta > Fasta.fasta
rm Fasta_AL_SR.fasta
##################################################################################################################################################
echo checking the sequnce lenght
bioawk -c fastx '{ print $name, length($seq) }' < Fasta.fasta |awk '{print$2}'|head -n 1 >fasta_lenght
cat fasta_lenght |while read LINE; do
if (($LINE % 3 == 0 )) 
then
{
##################################################################################################################################################
echo Infering approximately-maximum-likelihood phylogenetic trees from alignments
raxmlHPC -s Fasta.fasta -n out -m GTRCAT -f a -x 123 -N autoMRE -p 456 -T 40
mv RAxML_bestTree.out Tree.tree
rm Fasta.fasta.reduced RAxML_bipartitionsBranchLabels.out RAxML_bipartitions.out RAxML_bootstrap.out RAxML_info.out
###################################################################################################################################################
echo Runnning codeml for all branch 
#codeml codeml_Model0.ctl
slimcodeml codeml_Model0.ctl
echo codeml Run completed 
cat modelzero|grep 'lnL\|tree length =\|kappa (ts/tv) =\|omega (dN/dS) =' modelzero|xargs -n20 >modelzero.txt
echo Step Two
echo Species from the input 
grep "^>" Fasta.fasta |sed 's/>//g' >Species
cat Species 
echo Running codeml for each branch 
cat Species |while read r;do sed '1s/^/outfile ='$r'\n/' codeml_Model2.ctl > Model2.ctl|sed -E 's/\<('"$r"':[0-9.]*)/\1#1/' Tree.tree > T.Tree|slimcodeml Model2.ctl;done
echo Codeml Run completed
cat Species|while read r;do grep 'lnL\|tree length =\|kappa (ts/tv) =\|w (dN/dS) for branches' $r|xargs -n20| paste -|paste - modelzero.txt |sed 's/ntime:/ntime: /g'|sed 's/np:/np: /g'|awk 'BEGIN {OFS="\t";print "Species", "lnL_model0", "omega_model0", "lnL_model2", "omega_background", "omega_FW", "LRT", "Trend"};{gsub("lnL:",""); gsub(".fasta","");if ($20>$19) {print "'$r'", $25, $38, $5, $19, $20, 2*($5-$25), "faster"} else {print "'$r'",$25, $38, $5, $19, $20, 2*($5-$25), "slower"}}';done |awk '!NF || !seen[$0]++' >Results.txt
#add awk '$(NF-1) >3.84' for significance genes
echo Outputting the results
mv Results.txt Results_$(echo "$1" | cut -f 1 -d '.').txt
mv Tree.tree Tree_$(echo "$1" | cut -f 1 -d '.').tree
rm 2NG.dN 2NG.dS 2NG.t 4fold.nuc lnf Model2.ctl modelzero modelzero.txt rst rst1 rub T.Tree Fasta.fasta
cat Species|while read r;do rm $r ;done
rm Species
rm fasta_lenght
}
##################################################################################################################################################
else
{
   echo "Your seq is not mode 3" >> Results_$(echo "$1" | cut -f 1 -d '.').txt
}
fi
done
