# Cuaderno de Proceso: Análisis de Datos de Secuenciación de Bisulfito con Bismark

Este repositorio contiene el flujo de trabajo y la documentación para el análisis de datos de secuenciación con bisulfito (BS-Seq / RRBS).

Debido a que durante la secuenciación por bisulfito las citosinas no metiladas se leen como **T**, mientras que las **C** metiladas permanecen como **C**, se requiere un alineador especializado como **Bismark**. Esta herramienta alinea las lecturas tratadas con bisulfito respecto a una secuencia de referencia para determinar el estado de metilación.

---

## 📁 Estructura del Proyecto y Directorios

El flujo de trabajo se organiza entre el disco interno (procesamiento activo) y un disco externo (almacenamiento de resultados):

```text
📁 Disco interno (Procesamiento)
├── 📁 /home/mariana_castillo/
│   ├── 📄 Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz  (Genoma de referencia)
│   ├── 📁 Bisulfite_Genome/                               (Índice de Bismark)
│   └── 📁 datos/
│       └── 📁 trimmed/                                    (Secuencias .fq.gz de entrada)
│
📁 Disco externo (Resultados)
└── 📁 /mnt/d/
    └── 📁 Resultados/
        ├── 📁 alineamiento/                               (Archivos .bam alineados)
        └── 📁 qc_alineamiento/                            (Reportes .txt de Bismark)
```

---

## 📊 Dataset de Estudio

Para este análisis se utilizó el conjunto de datos del artículo:
> **Huang et al. (2021)**. *Effect of concentration and duration of particulate matter exposure on the transcriptome and DNA methylome of bronchial epithelial cells*. GEO: [GSE155615](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155615).

Según lo registrado en **SRA Run Selector**, las muestras provienen de secuenciación **RRBS** con lectura *Single-end*.

---

## 🚀 Flujo de Trabajo

### 1. Descargar los datasets de GEO

```bash
prefetch SRR12375071
```

### 2. Conversión de `.sra` a `.fastq`

Los archivos descargados en formato `.sra` deben convertirse a `.fastq`, que es el formato requerido por las herramientas posteriores.

```bash
# 1. Definir y crear la carpeta de destino
DESTINO="../datos/datosfastq"
mkdir -p "$DESTINO"
echo "Los archivos FASTQ se guardarán en: $DESTINO"
echo "-------------------------------------------"

find . -name "*.sra" -type f | while read archivo; do
    base=$(basename "$archivo" .sra)
    
    # Verificar existencia del archivo
    if [ ! -f "$archivo" ]; then
        echo "ERROR: No existe $archivo"
        continue
    fi
    
    # Verificar integridad
    echo "Verificando integridad de $base..."
    vdb-validate "$archivo"
    
    if [ $? -eq 0 ]; then
        echo "Convirtiendo: $base"
        date
        
        # 2. Convertir y enviar archivos directamente a DESTINO
        fasterq-dump -e 4 "$archivo" -O "$DESTINO/"
        
        # 3. Comprimir los archivos generados
        if [ -f "$DESTINO/${base}_1.fastq" ] || [ -f "$DESTINO/${base}.fastq" ]; then
            pigz -p 10 "$DESTINO/${base}"*.fastq
            echo "✓ $base convertido y comprimido en $DESTINO"
        else
            echo "✗ Error en la conversión de $base (no se encontraron archivos fastq en $DESTINO)"
        fi
    else
        echo "✗ $base está corrupto, saltando..."
    fi
    
    echo "---"
done
```

### 3. Control de Calidad Previo (FastQC)

Se evalúa la calidad de los datasets descargados antes del alineamiento:

```bash
fastqc datos/datosfastq/*.fastq.gz -o datos/QC_PRE/
```

### 4. Trimming y Optimización de Lecturas

**Trim Galore** elimina adaptadores y bases de baja calidad. Cuenta con un modo específico para RRBS que elimina artefactos de reparación de extremos y sesgos de metilación en librerías digeridas con *MspI*.

Adicionalmente, se utiliza **Clumpify** para ordenar los grupos de reads por similitud, lo cual reduce sensiblemente el tamaño de los archivos comprimidos finales y mejora la eficiencia.

```bash
# Trimming con Trim Galore (modo RRBS)
trim_galore --rrbs --fastqc --cores 10 -o ~/Resultados/trimmed/ ./datos/datosfastq/*.fastq.gz

# Reducción de espacio y duplicados con Clumpify
for file in ~/Resultados/trimmed/*_trimmed.fq.gz; do
    nombre_base=$(basename "$file")
    clumpify.sh in="$file" out="~/Resultados/clumpified/$nombre_base" dedupe=t subs=2
done
```

---

## 🧬 Alineamiento con Bismark

### Instalación de Bismark

Se puede instalar utilizando `mamba` o `micromamba`:

```bash
mamba install -c bioconda bismark
# O con micromamba:
micromamba activate
```

### Generación del Índice

Requiere como parámetro el directorio donde se encuentra el genoma de referencia en formato FASTA:

```bash
bismark prepare --parallel 6 .
```

### Configuración del Ambiente y Alineamiento

Creación del entorno específico de ejecución:

```bash
micromamba create -n bismark_env -c conda-forge -c bioconda bismark bowtie2 perl
```

Sintaxis general:

```bash
# Paired-end
bismark [options] --genome <genome_folder> {-1 <mates1> -2 <mates2> | <singles>}

# Single-end
bismark --genome . -p 3 -1 R1.fq.gz 

# Solución de ruta PERL5LIB si ocurren errores de compatibilidad
PERL5LIB="$CONDA_PREFIX/lib/perl5/core_perl:$CONDA_PREFIX/lib/perl5/5.32/core_perl:$CONDA_PREFIX/lib/perl5/5.32/site_perl:$CONDA_PREFIX/lib/perl5/site_perl:$CONDA_PREFIX/lib/perl5/5.32/vendor_perl:$CONDA_PREFIX/lib/perl5/vendor_perl" bismark --genome . -p 3 SRR12375071/SRR12375071.fastq.gz
```

Bucle principal de alineamiento y transferencia de resultados:

```bash
for archivo in ~/datos/trimmed/*.fq.gz; do
    echo "========================================"
    echo "Alineando secuencia: $archivo"
    
    PERL5LIB="$CONDA_PREFIX/lib/perl5/core_perl:$CONDA_PREFIX/lib/perl5/5.32/core_perl:$CONDA_PREFIX/lib/perl5/5.32/site_perl:$CONDA_PREFIX/lib/perl5/site_perl:$CONDA_PREFIX/lib/perl5/5.32/vendor_perl:$CONDA_PREFIX/lib/perl5/vendor_perl" bismark --genome . -p 3 -o Resultados/alineamiento/ "$archivo"
    
    if [ $? -eq 0 ]; then
        echo "Alineamiento exitoso."
        echo "Enviando .bam a disco externo..."
        mv Resultados/alineamiento/*.bam /mnt/d/Resultados/alineamiento
        echo "Enviando reporte .txt a qc_alineamiento..."
        mv /mnt/d/Resultados/alineamiento/*_report.txt /mnt/d/Resultados/qc_alineamiento/
        echo "Eliminando $archivo original para liberar espacio..."
        rm "$archivo"
    else
        echo "ALERTA: Error procesando $archivo. Conservando archivos por seguridad."
    fi
done

echo "========================================"
echo "¡Proceso finalizado!"
```

---

## 📈 Evaluación de Resultados de Alineamiento

**Eficiencia obtenida:** `67.9%`

> **Nota sobre el rendimiento:** Se obtuvo una eficiencia de alineamiento de 67.9%. Debido a que el ADN sufre daños químicos durante la conversión con bisulfito, es esperable que la amplificación y el alineamiento sean más difíciles (valores típicos aceptables oscilan entre 70% y 80%).

### Pruebas de Diagnóstico

Para evaluar si la baja eficiencia se debía a problemas de direccionalidad o al tipo de lectura (*paired-end* vs *single-end*), se realizó un muestreo de lecturas.

#### 1. Muestreo con SeqKit

Se extrajo el 10% del dataset original para realizar pruebas rápidas:

```bash
seqkit sample -p 0.1 input.fq.gz -o subset.fq.gz
```

#### 2. Alineamiento Direccional vs No Direccional

```bash
# Alineamiento no direccional
PERL5LIB="$CONDA_PREFIX/lib/perl5/core_perl:$CONDA_PREFIX/lib/perl5/5.32/core_perl:$CONDA_PREFIX/lib/perl5/5.32/site_perl:$CONDA_PREFIX/lib/perl5/site_perl:$CONDA_PREFIX/lib/perl5/5.32/vendor_perl:$CONDA_PREFIX/lib/perl5/vendor_perl" bismark --genome ./ref/ --non_directional -p 3 SRR12375071/SRR12375071.fastq.gz
```

**Conclusión:** Tanto el modo direccional como el no direccional arrojaron exactamente la misma eficiencia (**67.9%**). Sin embargo, el modo direccional requirió un tiempo de cómputo superior.

---

## 🔍 Extracción de Metilación y QC

### Control de Calidad: M-bias

Para detectar posibles sesgos de metilación según la posición de la base dentro de la lectura:

```bash
for bam in /mnt/d/Resultados/alineamiento/*.bam; do
    echo "Extrayendo M-bias de: $bam"
    bismark_methylation_extractor -s --mbias_only --multicore 3 -o /mnt/d/Resultados/qc_alineamiento/ "$bam"
done
```

### Reporte Consolidado con MultiQC

```bash
# Creación e instalación en entorno independiente
micromamba create -n multiqc_env python=3.10 -c conda-forge -y
pip install multiqc

# Generación del reporte
multiqc /mnt/d/Resultados/qc_alineamiento/ -o multiqc/
```

### Extracción de Contextos de Metilación (`bismark extract`)

La extracción discrimina de forma independiente entre citosinas en contexto **CpG**, **CHG** y **CHH**. En las lecturas resultantes se incluyen los identificadores de cadena origen:
* `OT`: Original top strand
* `CTOT`: Complementary to original top strand
* `OB`: Original bottom strand
* `CTOB`: Complementary to original bottom strand

#### Formato bedGraph
Genera salidas organizadas por coordenadas genómicas (inicio basado en 0, fin basado en 1) específicas para citosinas en contexto CpG.

```bash
for bam in /mnt/d/Resultados/alineamiento/*.bam; do
    echo "Extrayendo metilación de: $bam"
    bismark extract \
        --ignore 1 \
        --gzip \
        --bedGraph \
        --counts \
        --buffer_size 10G \
        --multicore 3 \
        -o Resultados/matriz/ \
        "$bam"
done
```

**Resultado:** Archivos de conteo `.cov.gz` (matriz de metilación) listos para importar a entornos de análisis posterior en **R** (p. ej., `methylKit` o `bsseq`).

---

## 📚 Bibliografía

* **Huang, S. K., Tripathi, P., Koneva, L. A., Cavalcante, R. G., Craig, N., Scruggs, A. M., Sartor, M. A., Deng, F., & Chen, Y. (2021).** *Effect of concentration and duration of particulate matter exposure on the transcriptome and DNA methylome of bronchial epithelial cells*. Environmental Epigenetics, 7(1), dvaa022. https://doi.org/10.1093/eep/dvaa022
* **Krueger, F. (2015).** *Trim Galore!: A wrapper around Cutadapt and FastQC to consistently apply adapter and quality trimming to FastQ files, with extra functionality for RRBS data*. Babraham Bioinformatics. https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/