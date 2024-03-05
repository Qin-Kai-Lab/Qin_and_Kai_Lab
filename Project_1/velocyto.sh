velocyto run10x  ./control /home/flyhunter/Wang/data/reference/refdata-gex-mm10-2020-A/genes/genes.gtf \
--samtools-threads 16 \
--samtools-memory 4000 


velocyto run10x ./control /home/flyhunter/Wang/data/reference/refdata-gex-mm10-2020-A/genes/genes.gtf --samtools-threads 12 --samtools-memory 4000

velocyto run10x ./stress /home/flyhunter/Wang/data/reference/refdata-gex-mm10-2020-A/genes/genes.gtf \
--samtools-threads 12 --samtools-memory 4000