# header of guideList:
# Guide_Name,Prefix,Guide_Sequence,Suffix,Target_Gene
# GAPDH_1,TTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCG,GCGGTGAGAAGCGCAGTCGG,GTTTAAGAGCTATGCTGGAAACAGCATAGCAAG,GAPDH

# for gRNAs of interest
awk -F '\t' '{if($4==3) print $1","toupper("tttcgatttcttggctttatatatcttgtggaaaggacgaaacaccG")","$2","toupper("GTTTaagagctatgctggaaacagcatagcaagt")","$1}' gRNAs_all.08132024_FINAL.txt >> guideList.csv

# for positive controls
awk 'BEGIN{FS=OFS="\t"}{if($3=="Positive_controls" && $NF==3) print $0}' gRNAs_all.08132024_FINAL.txt |sort -k 1 > gRNAs_positive.txt
wc -l gRNAs_positive.txt # 66/3=22

touch order.txt
for i in {1..22}
do
	echo "1
2
3" >> order.txt
done

paste gRNAs_positive.txt order.txt |awk -F '\t' '{print $1"_"$NF","toupper("tttcgatttcttggctttatatatcttgtggaaaggacgaaacaccG")","$2","toupper("GTTTaagagctatgctggaaacagcatagcaagt")","$1}' >> guideList.csv

rm order.txt

# for negative controls
awk 'BEGIN{FS=OFS="\t"}{if($3=="Negative_controls" && $NF==3) print $0}' gRNAs_all.08132024_FINAL.txt |sort -k 1 > gRNAs_negative.txt
wc -l gRNAs_negative.txt #60

touch order.txt
for i in {1..60}
do
        echo "${i}" >> order.txt
done

paste gRNAs_negative.txt order.txt |awk -F '\t' '{print $1","toupper("tttcgatttcttggctttatatatcttgtggaaaggacgaaacaccG")","$2","toupper("GTTTaagagctatgctggaaacagcatagcaagt")",ntc_"$NF}' |sed -rn 's/NonTargetingControlGuideForHuman/ntc/g' >> guideList.csv  

########## extract unique guides
awk -F ',' 'NR>1{print $3}' guideList.csv |sort |uniq -c |awk '{if($1==1)print $2}' > uniq_guides.list
awk -F ',' 'NR>1{print $3}' guideList.csv |sort |uniq -c |awk '{if($1!=1)print $2}' > dup_guides.list

head -n 1 guideList.csv > guideList.uniq.csv
touch tmp.sh tmp.txt
for i in `cat uniq_guides.list`
do
	echo "awk 'BEGIN{FS=OFS=\",\"}{if(\$3==\"${i}\")print \$0}' guideList.csv >> guideList.uniq.csv" >> tmp.sh
done
for i in `cat dup_guides.list`
do 
	echo "awk 'BEGIN{FS=OFS=\",\"}{if(\$3==\"${i}\")print \$0}' guideList.csv > tmp.txt
tail -n 1 tmp.txt >> guideList.uniq.csv
rm tmp.txt" >> tmp.sh
done

sh tmp.sh
rm tmp.sh tmp.txt

# add chromosome info
awk -F '\t' 'NR==FNR{a=$3;b[a]=$0;next}{OFS="\t";c=$28;d[c]=$8;if(b[c]){print b[c],d[c]}}' guideList.uniq.txt /scratch/hx37930/project/MTAG/06.crispr/results/gRNAs/final_gRNAs/positive_control_22genes.gRNAs.tsv > /scratch/hx37930/project/MTAG/06.crispr/sceptre/yk3_19117_35/guideList.chr.txt
awk -F '\t' 'NR==FNR{a=$3;b[a]=$0;next}{OFS="\t";c=$26;d[c]=$5;if(b[c]){print b[c],d[c]}}' guideList.uniq.txt /scratch/hx37930/project/MTAG/06.crispr/results/gRNAs/final_gRNAs/fine_mappedSNP.PIP05.removeTwoRegions.gRNAs.tsv >> /scratch/hx37930/project/MTAG/06.crispr/sceptre/yk3_19117_35/guideList.chr.txt
awk -F '\t' '{if($1~/ntc_/)print $0"\tNA:NA-NA"}' guideList.uniq.csv |sed 's/,/\t/g' >> /scratch/hx37930/project/MTAG/06.crispr/sceptre/yk3_19117_35/guideList.chr.txt
# remove duplicate target
awk -F '\t' 'NR==FNR{a=$1;b[a]=$0;next}{OFS="\t";c=$1;d[c]=$0;if(b[c]){print b[c]}}' /scratch/hx37930/project/MTAG/06.crispr/sceptre/yk3_19117_35/guideList.chr.txt  guideList.uniq.txt > /scratch/hx37930/project/MTAG/06.crispr/sceptre/yk3_19117_35/guideList.chr.uniq.txt
sed -i 's/:/\t/g;s/-/\t/g' /scratch/hx37930/project/MTAG/06.crispr/sceptre/yk3_19117_35/guideList.chr.uniq.txt
sed -i '1i Guide_Name\tPrefix\tGuide_Sequence\tSuffix\tTarget_Gene\tchr' /scratch/hx37930/project/MTAG/06.crispr/sceptre/yk3_19117_35/guideList.chr.uniq.txt
