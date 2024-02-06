## Pipeline Genotipagem:
***Este pipeline nos permitirá investigar a genotipagem do vírus SARS-CoV-2. Ele foi fornecido durante o curso de pós graduação de Bioinformática aplicada a genômica médica do Instituto de Ensino e Pesquisa Albert Einstein, pelo Prof. Deyvid Amgarten.*** Aconselha-se executá-lo no Google Colaboratory.

- Remove a pasta sample_data porque ela vem por padrão do Google Colabs e não é necessária
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

- ORF1ab:YP_009724389.1:exon1:c.C2772T:**p.F924F**,
- ORF1a:YP_009725295.1:exon1:c.C2772T:**p.F924F**,
- nsp3:YP_009742610.1:exon1:c.C318T:**p.F106F**
- ORF1ab:YP_009724389.1:exon1:c.C5247T:**p.N1749N**,
- ORF1a:YP_009725295.1:exon1:c.C5247T:**p.N1749N**,
- nsp3:YP_009742610.1:exon1:c.C2793T:**p.N931N**,
- ORF1ab:YP_009724389.1:exon1:c.C5775T:**p.F1925F**,
- ORF1a:YP_009725295.1:exon1:c.C5775T:**p.F1925F**,
- nsp3:YP_009742610.1:exon1:c.C3321T:**p.F1107F**,
- ORF1ab:YP_009724389.1:exon1:c.G7671T:**p.A2557A**,
- ORF1a:YP_009725295.1:exon1:c.G7671T:**p.A2557A**,
- nsp3:YP_009742610.1:exon1:c.G5217T:**p.A1739A**,
- ORF1ab:YP_009724389.1:exon1:c.T9378A:**p.I3126I**,
- ORF1a:YP_009725295.1:exon1:c.T9378A:**p.I3126I**,
- nsp4:YP_009725300.1:exon1:c.T1089A:**p.I363I**,
- ORF1ab:YP_009724389.1:exon2:c.A13448G:**p.K4483R**,
- nsp12:YP_009725307.1:exon2:c.A272G:**p.K91R**,
- ORF1ab:YP_009724389.1:exon2:c.G14026A:**p.D4676N**,
- nsp12:YP_009725307.1:exon2:c.G850A:**p.D284N**,
- ORF1ab:YP_009724389.1:exon2:c.C14144T:**p.P4715L**,
- nsp12:YP_009725307.1:exon2:c.C968T:**p.P323L**,
- S:YP_009724390.1:exon1:c.C67A:**p.Q23K**,
- S:YP_009724390.1:exon1:c.G769A:**p.G257S**,
- S:YP_009724390.1:exon1:c.A1841G:**p.D614G**,
- ORF3a:YP_009724391.1:exon1:c.C512T:**p.S171L**,
- ORF6:YP_009724394.1:exon1:c.96_98delinsTAC:**p.I33T**,
- N:YP_009724397.2:exon1:c.608_610delinsAAC:**p.R203_G204delinsKR**,
- N:YP_009724397.2:exon1:c.T875C:**p.I292T**.
![Screenshot from 2024-02-04 19-52-37](https://github.com/KairaCristina/TCC_metagen-mica_Instituto_Albert_Einstein/assets/131777938/98eb2d84-dc1f-422e-ae99-eeb62920f856)
Em destaque a posição 614 com o aminoácido glicina, PDB 7m8k.

https://github.com/KairaCristina/TCC_metagen-mica_Instituto_Albert_Einstein/assets/131777938/225cce7f-5e23-4eac-aaa1-2a1c81fbcc8a

PDB 7cn4 com "mutações" virtuais para os aminoácidos encontrados na amostra analisada acima.

