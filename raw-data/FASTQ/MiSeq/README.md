# Create soft links

Doing this makes it easier to use array jobs, like at `code/01_spaceranger/spaceranger.sh`. This is similar to https://github.com/LieberInstitute/DLPFC_snRNAseq/tree/main/raw-data. 

## Create dirs

Create directories with soft links to the fastq files:

```bash
while IFS= read -r SAMPLE
do
  echo "${SAMPLE}"
  mkdir ${SAMPLE}
done < "../../../code/01_spaceranger/samples_miseq.txt"
```

## Create soft links

Create the soft links based on the sample information from `raw-data/sample_info/Visium_HPC_Round1+2_110321_Master_SCP.xlsx`.

```bash
ln -s /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/1v-h* V10B01-085_A1/
ln -s /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/2v-h* V10B01-085_B1/
ln -s /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/3v-h* V10B01-085_C1/
ln -s /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/4v-h* V10B01-085_D1/
ln -s /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/5v-scp* V10B01-086_A1/
ln -s /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/v-scp_S6* V10B01-086_B1/
ln -s /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/7v-scp* V10B01-086_C1/
ln -s /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/8v-scp* V10B01-086_D1/
```

## Result

```bash
 $ ls -lh */
V10B01-085_A1/:
total 1.0K
lrwxrwxrwx 1 lcollado lieber_lcolladotor 121 Nov 18 12:19 1v-h_S1_L001_R1_001.fastq.gz -> /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/1v-h_S1_L001_R1_001.fastq.gz
lrwxrwxrwx 1 lcollado lieber_lcolladotor 121 Nov 18 12:19 1v-h_S1_L001_R2_001.fastq.gz -> /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/1v-h_S1_L001_R2_001.fastq.gz

V10B01-085_B1/:
total 1.0K
lrwxrwxrwx 1 lcollado lieber_lcolladotor 121 Nov 18 12:19 2v-h_S2_L001_R1_001.fastq.gz -> /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/2v-h_S2_L001_R1_001.fastq.gz
lrwxrwxrwx 1 lcollado lieber_lcolladotor 121 Nov 18 12:19 2v-h_S2_L001_R2_001.fastq.gz -> /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/2v-h_S2_L001_R2_001.fastq.gz

V10B01-085_C1/:
total 1.0K
lrwxrwxrwx 1 lcollado lieber_lcolladotor 121 Nov 18 12:19 3v-h_S3_L001_R1_001.fastq.gz -> /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/3v-h_S3_L001_R1_001.fastq.gz
lrwxrwxrwx 1 lcollado lieber_lcolladotor 121 Nov 18 12:19 3v-h_S3_L001_R2_001.fastq.gz -> /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/3v-h_S3_L001_R2_001.fastq.gz

V10B01-085_D1/:
total 1.0K
lrwxrwxrwx 1 lcollado lieber_lcolladotor 121 Nov 18 12:19 4v-h_S4_L001_R1_001.fastq.gz -> /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/4v-h_S4_L001_R1_001.fastq.gz
lrwxrwxrwx 1 lcollado lieber_lcolladotor 121 Nov 18 12:19 4v-h_S4_L001_R2_001.fastq.gz -> /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/4v-h_S4_L001_R2_001.fastq.gz

V10B01-086_A1/:
total 1.0K
lrwxrwxrwx 1 lcollado lieber_lcolladotor 123 Nov 18 12:19 5v-scp_S5_L001_R1_001.fastq.gz -> /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/5v-scp_S5_L001_R1_001.fastq.gz
lrwxrwxrwx 1 lcollado lieber_lcolladotor 123 Nov 18 12:19 5v-scp_S5_L001_R2_001.fastq.gz -> /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/5v-scp_S5_L001_R2_001.fastq.gz

V10B01-086_B1/:
total 1.0K
lrwxrwxrwx 1 lcollado lieber_lcolladotor 122 Nov 18 12:20 v-scp_S6_L001_R1_001.fastq.gz -> /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/v-scp_S6_L001_R1_001.fastq.gz
lrwxrwxrwx 1 lcollado lieber_lcolladotor 122 Nov 18 12:20 v-scp_S6_L001_R2_001.fastq.gz -> /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/v-scp_S6_L001_R2_001.fastq.gz

V10B01-086_C1/:
total 1.0K
lrwxrwxrwx 1 lcollado lieber_lcolladotor 123 Nov 18 12:19 7v-scp_S7_L001_R1_001.fastq.gz -> /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/7v-scp_S7_L001_R1_001.fastq.gz
lrwxrwxrwx 1 lcollado lieber_lcolladotor 123 Nov 18 12:19 7v-scp_S7_L001_R2_001.fastq.gz -> /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/7v-scp_S7_L001_R2_001.fastq.gz

V10B01-086_D1/:
total 1.0K
lrwxrwxrwx 1 lcollado lieber_lcolladotor 123 Nov 18 12:19 8v-scp_S8_L001_R1_001.fastq.gz -> /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/8v-scp_S8_L001_R1_001.fastq.gz
lrwxrwxrwx 1 lcollado lieber_lcolladotor 123 Nov 18 12:19 8v-scp_S8_L001_R2_001.fastq.gz -> /dcs04/lieber/lcolladotor/libdRNASEQcore_LIBD001/raw-data/2021-11-17_MiSeq_HPC/bcl2fastq_out/8v-scp_S8_L001_R2_001.fastq.gz
```