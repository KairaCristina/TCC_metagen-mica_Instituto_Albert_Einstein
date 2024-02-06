## Pipeline Metagenômica:
***Este pipeline nos permitirá investigar vários micro-organismos presentes na amostra analisada. Ele foi fornecido durante o curso de pós-graduação em Bioinformática aplicada à genômica médica pelo Instituto de Ensino e Pesquisa Albert Einstein, ministrado pelo Prof. Deyvid Amgarten.*** Recomenda-se executá-lo no Google Colaboratory.

## **Parte 1** - Organização do espaço de trabalho
- Remove a pasta sample_data porque ela vem por padrão do Google Colab e não é necessária para este pipeline.

```
! rm -rf sample_data/
```
- Conexão do Google Drive: Será nosso repositório persistente de dados para não ter que recomeçar do zero as anotações.

```
from google.colab import drive
drive.mount('/content/drive')
```
## **Parte 2** - Configurar o ambiente para execução do pipeline
- Nesta etapa, vamos criar um ambiente conda e instalar as ferramentas que serão utilizadas no pipeline.
```
%%bash
# Instala o gerenciador de envs Conda
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b
```

```
import os
os.environ['PATH'] += ':/root/miniforge3/bin/'
```
- Instala todos os programas necessários
```
!mamba install --channel bioconda \
  fastqc cutadapt kraken2 krona bwa \
  samtools spades --yes
```
- testa se os programas foram instalados corretamente
```
!samtools
```
- Faz download dos arquivos de taxid necessários no krona tools
```
!ktUpdateTaxonomy.sh
```
## **Parte 3** - Download dos dados
- Download da amostra:
```
!wget https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_METAGENOMICA_S9_R1_001.fastq.gz
```
```
!wget https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_METAGENOMICA_S9_R2_001.fastq.gz
```
- Download kraken2 database
```
!wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20231009.tar.gz
```
- Download genoma humano (hospedeiro) para remoção de contaminantes
```
!wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
!tar -xf Homo_sapiens_UCSC_hg38.tar.gz
```
```
HOST = '/content/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'
```
- Reorganizar os arquivos em uma estrutura intelegível de pastas
```
!mkdir -p kraken2 kraken2-db fastq fastqc cutadapt host bwa
!tar -xf k2_standard_08gb_20231009.tar.gz --directory kraken2-db/
!rm k2_standard_08gb_20231009.tar.gz
!rm Miniforge3-Linux-x86_64.sh
!mv /content/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa host/
HOST = 'genome.fa'
!mv /content/patient_joao_METAGENOMICA_S9_R1_001.fastq.gz /content/fastq
!mv /content/patient_joao_METAGENOMICA_S9_R2_001.fastq.gz /content/fastq
```
- Coletar nome da amostra em uma variável
```
temp = !basename fastq/*_R1_* _R1_001.fastq.gz
SAMPLE = temp[0]
!echo {SAMPLE} {HOST}
```
## **Parte 4** - Controle de qualidade e limpeza das sequências
- Nesta etapa, vamos gerar um report de qualidade para os dados de sequenciamento da amostra e com base nele, vamos proceder com algumas limpezas (filtragem e trimagem).
- Gerar relatórios de qualidade do sequenciamento com fastqc
```
!fastqc fastq/{SAMPLE}_R1_001.fastq.gz \
  fastq/{SAMPLE}_R2_001.fastq.gz \
  -o fastqc/
```
### Imagem tabela de qualidade arquivo <patient_joao_METAGENOMICA_S9_R1_001_fastqc.html>
![Screenshot from 2024-02-06 16-27-19](https://github.com/KairaCristina/TCC_metagen-mica_Instituto_Albert_Einstein/assets/131777938/c21a9f28-5cf0-4788-a5b6-60da80695bd9)

### Imagem gráfico de qualidade arquivo <patient_joao_METAGENOMICA_S9_R1_001_fastqc.html>
![Screenshot from 2024-02-06 16-27-34](https://github.com/KairaCristina/TCC_metagen-mica_Instituto_Albert_Einstein/assets/131777938/cf03c8b0-0995-4842-beaa-ea0e0a64bce9)

### Imagem tabela de qualidade arquivo <patient_joao_METAGENOMICA_S9_R2_001_fastqc.html>
![Screenshot from 2024-02-06 16-26-29](https://github.com/KairaCristina/TCC_metagen-mica_Instituto_Albert_Einstein/assets/131777938/fed95fea-d469-4364-9443-06f6d54f0c55)

### Imagem gráfico de qualidade arquivo <patient_joao_METAGENOMICA_S9_R2_001_fastqc.html>
![Screenshot from 2024-02-06 16-26-49](https://github.com/KairaCristina/TCC_metagen-mica_Instituto_Albert_Einstein/assets/131777938/7d6ef36e-bc4b-4dfe-8993-6a68c06feba1)

- Filtragem e trimagem das sequências com cutadapt
```
!cutadapt -u 5 -U 5 -u -9 -U -9 -m 50  \
  -o cutadapt/{SAMPLE}_cleaned_R1.fastq.gz \
  -p cutadapt/{SAMPLE}_cleaned_R2.fastq.gz \
  fastq/{SAMPLE}_R1_001.fastq.gz \
  fastq/{SAMPLE}_R2_001.fastq.gz > cutadapt/summary_cutadapt.txt
```
### Imagem de relátório arquivo <summary_cutadapt.txt>
![Screenshot from 2024-02-06 16-25-50](https://github.com/KairaCristina/TCC_metagen-mica_Instituto_Albert_Einstein/assets/131777938/8bdca848-023e-42ea-b708-289b70809c45)
## **Parte 4** - Remover contaminantes do hospedeiro
- Nesta etapa, precisamos pegar as sequências limpas na etapa anterior e proceder com o mapeamento no genoma de referência humano. Dessa forma, poderemos remover estas sequências que mapeias contra o genoma humano e ficar com apenas o que interessa para as demais etapas. 
- Indexar a referência
```
!bwa index host/{HOST}
```
- Proceder com o alinhamento contra o genoma completo humano
```
!bwa mem host/{HOST} \
  cutadapt/{SAMPLE}_cleaned_R1.fastq.gz \
  cutadapt/{SAMPLE}_cleaned_R2.fastq.gz \
  | samtools view -b > bwa/{SAMPLE}_mapped_host.bam
```
- Gerar BAM de não mapeados
```
!samtools view -u -f 12 \
  -b bwa/{SAMPLE}_mapped_host.bam \
  | samtools sort -n > bwa/{SAMPLE}_unmapped_host.bam
```
- Gerar relatório com os counts de reads mapeados em humano
```
!samtools flagstat bwa/{SAMPLE}_mapped_host.bam \
  > bwa/{SAMPLE}_mapped_host_flagstat.txt
```
- Gerar FASTQ de reads não mapeados
```
!samtools fastq bwa/{SAMPLE}_unmapped_host.bam \
  -1 bwa/{SAMPLE}_unmapped_host_R1.fastq \
  -2 bwa/{SAMPLE}_unmapped_host_R2.fastq
```
```
!ls -lh bwa/
```
## **Parte 5** - Identificação taxonômica
- Nesta etapa, vamos pegar os reads filtrados e não mapeados no hospedeiro para proceder com a identificação taxonômica dos reads. Esta etapa faz uma busca dos reads contra um banco de dados de organismos conhecidos, utilizando heurísticas específicas da ferramenta Kraken2.
O kraken2 possui vários bancos disponíveis, com diferentes composições de genomas de organismos. Neste exemplo prático, escolhemos o banco Standard 8GB que é uma versão com Bacterias, Vírus, Archaea e Humano do NCBI Refseq. Os bancos do kraken2 geralmente utilizam muita memória RAM e isso pode ser uma limitação de execução em computadores mais simples.
```
!kraken2 -db kraken2-db/ \
	--report kraken2/{SAMPLE}_kraken_report.txt \
  --output kraken2/{SAMPLE}_kraken2_NT.out \
	--minimum-base-quality 20 \
  --paired bwa/{SAMPLE}_unmapped_host_R1.fastq \
  bwa/{SAMPLE}_unmapped_host_R2.fastq

```
- Gerar relatórios
```
!ktImportTaxonomy kraken2/{SAMPLE}_kraken2_NT.out \
  -q 2 -t 3 \
  -o kraken2/{SAMPLE}_classification.html
```
## **Parte final** - Inspecionar os resultados
- Nesta última etapa, cabe ao analista ou médico patologista/infectologista inspecionar os resultados gerados pelos relatórios de diversidade para identificar possíveis patógenos na amostra.
É importante ressaltar que o resultado da identificação taxonômica gerado pelo kraken2 não deve ser levado em consideração unicamente. O processo de laudamente deve envolver uma etapa posterior de validação do achado, que poderia remover falso-positivos gerados pelo Kraken2 ou mesmo encontrar patógenos que não foram apontado pela identificação do Kraken2.

### Imagem de gráfico resultado final da análise, arquivo <patient_joao_METAGENOMICA_S9_classification.html>
![Screenshot from 2024-02-06 16-30-41](https://github.com/KairaCristina/TCC_metagen-mica_Instituto_Albert_Einstein/assets/131777938/5566b1fd-6a13-4424-b1f2-685c0e5e9037)