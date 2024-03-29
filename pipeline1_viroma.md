## Pipeline Viroma:
***Este pipeline foi fornecido durante o curso de pós graduação de Bioinformática aplicada a genômica médica do Instituto de Ensino e Pesquisa Albert Einstein, pelo Prof. Deyvid Amgarten.***

### **Parte 1** - Organização do espaço de trabalho
- Remove a pasta sample_data porque ela vem por padrão do Google Colab e não é necessária para este pipeline.

```
! rm -rf sample_data/
```
```
PACIENTE_S21_R1 = "https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_VIROMA_S21_R1_001.fastq.gz"
PACIENTE_S21_R2 = "https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_VIROMA_S21_R2_001.fastq.gz"
```
### **Parte 2** - Instalando ferramentas e arquivos a serem utilizados
- Configuração do ambiente de execução: Nesta etapa, vamos criar um ambiente conda e instalar as ferramentas que serão utilizadas no pipeline.
```
%%bash
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
- Testa se os programas foram instalados corretamente
```
!samtools
```
- Faz download dos arquivos de taxid necessários no krona tools

```
!ktUpdateTaxonomy.sh
```
- Nesta etapa, nós iremos fazer download dos arquivos fastq.gz da amostra para dentro do nosso ambiente de execução. Também iremos fazer download do banco utilizado pelo kraken2 para identificação taxonômica de organismos.
- Clone gdown - faz download de links compartilhados do Google Drive
- Download R1

```
%%bash
git clone https://github.com/circulosmeos/gdown.pl.git
./gdown.pl/gdown.pl https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_VIROMA_S21_R1_001.fastq.gz PACIENTE1_VIROMA_S21_R1_001.fastq.gz
```
- Download R2
```
%%bash
git clone https://github.com/circulosmeos/gdown.pl.git
./gdown.pl/gdown.pl https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_VIROMA_S21_R2_001.fastq.gz PACIENTE1_VIROMA_S21_R2_001.fastq.gz
```
- Download kraken2 database
```
!wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20231009.tar.gz
```
- Download genoma do hospedeiro para remoção de contaminantes.

```
!wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
!tar -xf Homo_sapiens_UCSC_hg38.tar.gz
HOST = '/content/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'
```
- Reorganizar os arquivos em uma estrutura intelegível de pastas.

```
!mkdir -p fastq kraken2 kraken2-db fastqc cutadapt host bwa
!mv *fastq.gz fastq/
!tar -xf k2_standard_08gb_20231009.tar.gz --directory kraken2-db/
!rm k2_standard_08gb_20231009.tar.gz
!rm Miniforge3-Linux-x86_64.sh
!mv {HOST} host/
```
- Coletar nome da amostra em uma variável
```
temp = !basename fastq/*_R1_* _R1_001.fastq.gz
SAMPLE = temp[0]
```
- Verificar se o nome que aparece abaixo é o do paciente escolhido.
```
!echo {SAMPLE} {HOST}
```
### **Parte 3** - Controle de qualidade e limpeza das sequências
- Nesta etapa, vamos gerar um report de qualidade para os dados de sequenciamento da amostra e com base nele, vamos proceder com algumas limpezas (filtragem e trimagem).
- Gerar relatórios de qualidade do sequenciamento com fastqc.
```
!fastqc fastq/{SAMPLE}_R1_001.fastq.gz \
  fastq/{SAMPLE}_R2_001.fastq.gz \
  -o fastqc/
```
# Imagem relatório fastqc.html (R1)
![Screenshot from 2024-01-22 08-45-33](https://github.com/KairaCristina/TCC_metagen-mica_Instituto_Albert_Einstein/assets/131777938/8a116dca-b72d-459f-9bea-454cb89f7233)
![Screenshot from 2024-01-22 08-47-54](https://github.com/KairaCristina/TCC_metagen-mica_Instituto_Albert_Einstein/assets/131777938/baf24883-2eb9-4b45-9597-59ecfc08f9df)
# Imagem relatório fastqc.html (R2)
![image](https://github.com/KairaCristina/TCC_metagen-mica_Instituto_Albert_Einstein/assets/131777938/7e808635-2813-4f41-a6bd-5794d2834a74)
![image](https://github.com/KairaCristina/TCC_metagen-mica_Instituto_Albert_Einstein/assets/131777938/3a07b98c-484e-4db6-9037-948043a542dc)


- Filtragem e trimagem das sequências com cutadapt

```
!cutadapt -u 5 -U 5 -u -9 -U -9 -m 50 \
  -o cutadapt/{SAMPLE}_cleaned_R1.fastq.gz \
  -p cutadapt/{SAMPLE}_cleaned_R2.fastq.gz \
  fastq/{SAMPLE}_R1_001.fastq.gz \
  fastq/{SAMPLE}_R2_001.fastq.gz > cutadapt/summary_cutadapt.txt
```
### Arquivo <summary_cutadapt.txt>
![image](https://github.com/KairaCristina/TCC_metagen-mica_Instituto_Albert_Einstein/assets/131777938/28e9faaf-94d0-4fcd-a3ab-63b7d50d0162)
```
!cutadapt -h
```
### **Parte 4** - Remover contaminantes do hospedeiro.
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
- Gerar relatório com os counts de reads mapeados em humano.
```
!samtools flagstat bwa/{SAMPLE}_mapped_host.bam \
  > bwa/{SAMPLE}_mapped_host_flagstat.txt
```
# Relatório .txt
![image](https://github.com/KairaCristina/TCC_metagen-mica_Instituto_Albert_Einstein/assets/131777938/4ab18463-4e36-427a-bbe7-6545e878dad8)
- Gerar FASTQ de reads não mapeados
```
!samtools fastq bwa/{SAMPLE}_unmapped_host.bam \
  -1 bwa/{SAMPLE}_unmapped_host_R1.fastq \
  -2 bwa/{SAMPLE}_unmapped_host_R2.fastq
```
- Conferir arquivos

```
!ls -lh bwa/
```
### **Parte 5** - Identificação taxonômica.
- Nesta etapa, vamos pegar os reads filtrados e não mapeados no hospedeiro para proceder com a identificação taxonômica dos reads. Esta etapa faz uma busca dos reads contra um banco de dados de organismos conhecidos, utilizando heurísticas específicas da ferramenta Kraken2.
O kraken2 possui vários bancos disponíveis, com diferentes composições de genomas de organismos. Neste exemplo prático, escolhemos o banco Standard 8GB que é uma versão com Bacterias, Vírus, Archaea e Humano do NCBI Refseq. Os bancos do kraken2 geralmente utilizam muita memória RAM e isso pode ser uma limitação de execução em computadores mais simples. Todavia, uma versão que utiliza o máximo de 8GB de RAM é disponibilizada para computadores pequenos e esta versão está sendo utilizada nesta prática.

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
# Imagem resultado final total
![image](https://github.com/KairaCristina/TCC_metagen-mica_Instituto_Albert_Einstein/assets/131777938/70a4204d-9061-4fdb-aeca-92019140a4da)
# Imagem resultado final viroma
![image](https://github.com/KairaCristina/TCC_metagen-mica_Instituto_Albert_Einstein/assets/131777938/a756601c-9681-41a0-a20f-2205c63b2841)
# Imagem resultado final Pseudomonas
![image](https://github.com/KairaCristina/TCC_metagen-mica_Instituto_Albert_Einstein/assets/131777938/69183d1c-6804-41a4-b3f6-b8903976a5cc)

### **Inspecionar resultados**
- Nesta última etapa, cabe ao analista ou médico patologista/infectologista inspecionar os resultados gerados pelos relatórios de diversidade para identificar possíveis patógenos na amostra.
É importante ressaltar que o resultado da identificação taxonômica gerado pelo kraken2 não deve ser levado em consideração unicamente. O processo de laudamente deve envolver uma etapa posterior de validação do achado, que poderia remover falso-positivos gerados pelo Kraken2 ou mesmo encontrar patógenos que não foram apontado pela identificação do Kraken2.

## Etapa adicional de validação: identificação taxônomica por contigs.fasta
Nesta etapa faremos o download do arquivo do paciente s21, com a etapa de montagem, realizada pela ferramenta Spades, dos reads em contigs (sequências maiores) e então realizar a mesma etapa de identificação taxônomica.
```
!mkdir Assembly
```
- Clone gdown - faz download de links compartilhados do Google Drive e download do arquivo pelo link
```
%%bash
git clone https://github.com/circulosmeos/gdown.pl.git
./gdown.pl/gdown.pl https://aulas-pos-hiae-public-data.s3.sa-east-1.amazonaws.com/TCC-metagenomica/patient_joao_METAGENOMICA_contigs.fasta PATIENT_joao_METAGENOMICA_contigs.fasta
```
```
!mv /content/PATIENT_joao_METAGENOMICA_contigs.fasta /content/Assembly
```
- Nova classificação para os contigs montados (sequencias mais longas)
```
!kraken2 -db kraken2-db/ \
	--report kraken2/{SAMPLE}_kraken_report_contigs.txt \
  --output kraken2/{SAMPLE}_kraken2_NT_contigs.out \
 /content/Assembly/PATIENT_joao_METAGENOMICA_contigs.fasta
```
- Gerar novo relatório, agora para os contigs
```
!ktImportTaxonomy kraken2/{SAMPLE}_kraken2_NT_contigs.out \
  -q 2 -t 3 \
  -o kraken2/patient_joao_METAGENOMICA_classification_contigs.html
```
- imagem arquivo <patient_joao_METAGENOMICA_classification_contigs.html>
![Screenshot from 2024-02-04 17-33-14](https://github.com/KairaCristina/TCC_metagen-mica_Instituto_Albert_Einstein/assets/131777938/a3e965c9-68de-4011-8cc5-bd475a2ec60d)

![Screenshot from 2024-02-04 17-33-42](https://github.com/KairaCristina/TCC_metagen-mica_Instituto_Albert_Einstein/assets/131777938/f034df5e-a2fc-4ca8-b177-51df611cfad8)

### Conclusão da etapa adicional:
A sugestão de montar os reads em contigs antes de realizar a identificação taxonômica pode oferecer várias vantagens em comparação com a busca direta dos reads. Aqui estão algumas das possíveis vantagens desse procedimento:

- Melhor resolução de sequências: Montar os reads em contigs pode resultar em sequências mais longas e contínuas, proporcionando uma melhor resolução do genoma. Isso é particularmente útil para identificar genes completos, regiões regulatórias e outras características genômicas que podem ser fragmentadas quando visualizadas apenas como reads individuais.

- Redução de artefatos de sequenciamento: A montagem em contigs pode ajudar a eliminar ou reduzir artefatos de sequenciamento, como erros de leitura e duplicatas, proporcionando uma representação mais precisa do genoma.

- Identificação de novos genes ou elementos genéticos: A montagem em contigs pode revelar novos genes ou elementos genéticos que podem não ser facilmente identificados quando apenas os reads são analisados. Isso é especialmente importante em estudos de genomas não completamente sequenciados ou em organismos pouco caracterizados.

