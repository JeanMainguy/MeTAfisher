# faa_file=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Acry_Genome/sequence.faa
# fna_file=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Acry_Genome/sequence.fna
# contig_file=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Acry_Genome/sequence.fasta
# faa_file_transformed=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Acry_Genome/sequence_transformed.faa
# fna_file_transformed=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Acry_Genome/sequence_transformed.fna
# faa_file=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Actino_pleuro/sequence.faa
# fna_file=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Actino_pleuro/sequence.fna
#
# faa_file_transformed=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Actino_pleuro/sequence_transformed.faa
# fna_file_transformed=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Actino_pleuro/sequence_transformed.fna
# contig_file=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Actino_pleuro/sequence.fasta
faa_file=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Alkalilimnicola_ehrlichii_MLHE-1/sequence.faa
fna_file=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Alkalilimnicola_ehrlichii_MLHE-1/sequence.fna

faa_file_transformed=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Alkalilimnicola_ehrlichii_MLHE-1/sequence_transformed.faa
fna_file_transformed=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Alkalilimnicola_ehrlichii_MLHE-1/sequence_transformed.fna
contig_file=/home/jean/Documents/Master_Bioinfo/S2_Bilbao/StageBAO/ProgamV2/data/Alkalilimnicola_ehrlichii_MLHE-1/sequence.fasta

contig=$(head -1 ${contig_file}|cut -d' ' -f1|cut -d'>' -f2)
contig=CP000453.1
echo $contig
cat ${fna_file} | sed -e "s/>lcl|${contig}[A-Za-z_]*/>${contig}|/g"|sed 's/_[0-9]*//'| sed 's/|0*/|/' > ${fna_file_transformed}
# cat tmp.txt

cat ${faa_file} | sed -e "s/>lcl|${contig}[A-Za-z_]*/>${contig}|/g"|sed 's/_[0-9]*//'| sed 's/|0*/|/' > ${faa_file_transformed}
# cat tmp.txt
