## Pipeline Genotipagem:
***Este pipeline foi fornecido durante o curso de pós graduação de Bioinformática aplicada a genômica médica do Instituto de Ensino e Pesquisa Albert Einstein, pelo Prof. Deyvid Amgarten.*** Aconselha-se executá-lo no Google Colaboratory.

- Remove a pasta sample_data porque ela vem porpadrão do Google Colabs e não é necessária
```
! rm -rf sample_data/
```
- Utilize o comando abaixo para instalar os programas necessários para executar o pipeline de genotipagem viral.
```
%%bash
sudo apt install tree fastqc bwa samtools bedtools freebayes
pip install cutadapt
```
- Download dos arquivos necessários da amostra (FASTQ)
```
PACIENTE_R1 = "https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_VIROMA_S21_R1_001.fastq.gz"
PACIENTE_R2 = "https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_VIROMA_S21_R2_001.fastq.gz"
```
```
%%bash
git clone https://github.com/circulosmeos/gdown.pl.git
./gdown.pl/gdown.pl https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_VIROMA_S21_R1_001.fastq.gz PACIENTE1_VIROMA_S21_R1_001.fastq.gz
```
```
%%bash
git clone https://github.com/circulosmeos/gdown.pl.git
./gdown.pl/gdown.pl https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_VIROMA_S21_R2_001.fastq.gz PACIENTE1_VIROMA_S21_R2_001.fastq.gz
```
- Download SARS-COV-2 Reference
```
%%bash
wget -nv https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/NC_045512.fasta
```
- Organizar arquivos nas suas respectivas pastas
```
%%bash
mkdir -p logs fastq passedQC reference mapped variants coverage
mv *fastq.gz fastq/
mv NC_045512.fasta reference/
```
- Coletar nome da amostra em uma variável
```
temp = !basename fastq/*_R1* _R1_001.fastq.gz
SAMPLE = temp[0]
print(SAMPLE)
```
- Controle de Qualidade das Sequências com FASTQC
```
!fastqc fastq/{SAMPLE}_R1_001.fastq.gz fastq/{SAMPLE}_R2_001.fastq.gz
```
- Limpeza das sequências com cutadapt
```
!cutadapt -u 5 -U 5 -u -9 -U -9 -m 50 -o passedQC/{SAMPLE}_cleaned_R1.fastq -p passedQC/{SAMPLE}_cleaned_R2.fastq fastq/{SAMPLE}_R1_001.fastq.gz fastq/{SAMPLE}_R2_001.fastq.gz
```
- Mapeamento dos reads no genoma de referência do SARS-CoV-2: Primeiramente, é necessário preparar a referência para o alinhamento. É o que chamamos de indexação da referência.
```
%%bash
bwa index reference/NC_045512.fasta
```
- Após a indexação, podemos usar o comando que faz o mapeamento chamado bwa mem. Lembrando sempre que o programa precisa dos fastqs limpos e da referência como entrada.
```
!bwa mem reference/NC_045512.fasta passedQC/{SAMPLE}_cleaned_R1.fastq passedQC/{SAMPLE}_cleaned_R2.fastq > mapped/{SAMPLE}_mapped_sarscov2.sam
```
- Organizar o SAM e gerar o BAM
```
!samtools sort mapped/{SAMPLE}_mapped_sarscov2.sam -o mapped/{SAMPLE}_mapped_sarscov2_sorted.bam
```
- Gerar dados de cobertura para o mapeamento
```
!bedtools bamtobed -i mapped/{SAMPLE}_mapped_sarscov2_sorted.bam > coverage/mapped_sarcov2.bed
!bedtools merge -i coverage/mapped_sarcov2.bed >coverage/mapped_sarscov2_merged.bed
!bedtools sort -i coverage/mapped_sarscov2_merged.bed >coverage/mapped_sarscov2_merged_sorted.bed
```
- Após gerar este .BED, vamos utilizar o comando bedtools coverage para gerar a cobertura média nesta região contínua.
```
!bedtools coverage -a coverage/mapped_sarscov2_merged_sorted.bed \
-b mapped/{SAMPLE}_mapped_sarscov2_sorted.bam -mean \
>coverage/results_coverage.bed
```
- Fazer a chamada de variantes com o Freebayes
```
!freebayes -f reference/NC_045512.fasta -p 1 mapped/{SAMPLE}_mapped_sarscov2_sorted.bam > variants/{SAMPLE}_variants.vcf
```
- Anotação do VCF com base no genoma de referência de SARS-CoV-2: Download do annovar.
```
!wget -nv https://github.com/Varstation/T1-2020/raw/master/annovar/annovar.zip
!wget -nv http://www.openbioinformatics.org/annovar/download/NC_045512v2_avGene.txt.gz
!wget -nv http://www.openbioinformatics.org/annovar/download/NC_045512v2_avGeneMrna.fa.gz
```
- Descompartar e remover os arquivos zipados
```
%%bash
unzip annovar.zip
gunzip -c NC_045512v2_avGene.txt.gz > reference/NC_045512v2_avGene.txt
gunzip -c NC_045512v2_avGeneMrna.fa.gz > reference/NC_045512v2_avGeneMrna.fa
rm annovar.zip *.gz
```
- Comando anotação
```
!sed -i 's/NC_045512.2/NC_045512v2/g' variants/{SAMPLE}_variants.vcf
```
```
!annovar/table_annovar.pl -buildver NC_045512v2 -vcfinput variants/{SAMPLE}_variants.vcf reference/ -protocol avGene -operation g --polish
```
- Filtrar as linhas com AF=1 para observar as variantes presentes no vírus SARS-CoV-2 presente nesta amostra.
```
!head -n 1 variants/{SAMPLE}_variants.vcf.NC_045512v2_multianno.txt > variants/{SAMPLE}_final_variants_annotated.txt
!grep ';AF=1;' variants/{SAMPLE}_variants.vcf.NC_045512v2_multianno.txt >> variants/{SAMPLE}_final_variants_annotated.txt
```
No arquivo <PACIENTE1_VIROMA_S21_final_variants_annotated.txt>, ao abri-lo no Excel ou LibreOffice Calc, você poderá visualizar, na coluna 'AAChange.avGene', a alteração de aminoácidos correspondente ao genoma do vírus presente na amostra analisada em cada anotação.

### Trocas de aminoácidos da amostra analisada:
AAChange

- F106F
- N931N
- F1107F
- A1739A
- I363I
- K91R
- D284N
- P323L
- Q23K
- G257S
- D614G
- S171L
- I33T
- G204delinsKR
- I292T
