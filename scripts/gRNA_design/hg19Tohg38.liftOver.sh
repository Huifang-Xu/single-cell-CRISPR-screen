prefix="fine_mappedSNP.all.removeTwoRegions"
reference="/scratch/hx37930/project/Hawaiian_usc/reference/hg19ToHg38.over.chain.gz"

cd /scratch/hx37930/project/MTAG/06.crispr/data 

/home/hx37930/liftOver ${prefix}.hg19.sort.txt ${reference} ${prefix}.hg19Tohg38.map.txt ${prefix}.hg19Tohg38.unmap.txt

# add header
sed -i '1i hg38chr\thg38start\thg38end\tSNP:hg19chr:hg19BP:A1:A2' ${prefix}.hg19Tohg38.map.txt

