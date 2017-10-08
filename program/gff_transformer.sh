# gff_file=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Acry_Genome/sequence.gff
# gff_file=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Acry_Genome/tmp_head.gff
# gff_file_transformed=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Acry_Genome/sequence_transformed.gff
gff_file=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Actino_pleuro/sequence.gff
gff_file_transformed=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Actino_pleuro/sequence_transformed.gff
gff_file=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Alkalilimnicola_ehrlichii_MLHE-1/sequence.gff
gff_file_transformed=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Alkalilimnicola_ehrlichii_MLHE-1/sequence_transformed.gff

rm ${gff_file_transformed}
while read p; do
    type=$(echo $p | cut -d$' ' -f3)

    # echo $type
    if [ ${type} == 'gene' ]
    then
        gene_name=$(echo $p |sed -n 's/.*ID=\(gene[0-9]*\).*/\1/p')
        nbg=$(echo $p | sed -n 's/.*locus_tag=[A-Za-z]*_\([0-9]*\).*/\1/p')

fi
    if [ ${type} == 'CDS' ]
    then
        contig=$(echo $p | cut -d$' ' -f1)
        nb=$(echo $p |sed -n 's/.*ID=cds\([0-9]*\).*/\1/p')
        nbi=$(echo $p |sed -n 's/.*Name=[A-Za-z]*_\([0-9]*\).*/\1/p')
        # echo $p
        echo $nbi
        # nbi=$(($nb+1))
        parent=$(echo $p | sed -n 's/.*Parent=\(gene[0-9]*\).*/\1/p')
        # if [ ${parent} == $gene_name]
        # then
            np=$(echo $p | sed -e "s/ID=cds${nb}/ID=${contig}|${nbg}/g")

            npt=$(echo -e $np |sed -e 's/ [ ]*/\t/g')
            echo -e $npt >> ${gff_file_transformed}
    # fi
fi
done <${gff_file}
cat ${gff_file_transformed} | sed -e 's/ [ ]*/\t/g'| sed 's/|0*/|/' > tmp.txt
mv tmp.txt ${gff_file_transformed}
