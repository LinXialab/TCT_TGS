import csv
import os
import pandas as pd
import csv
from multiprocessing import Pool

############################################
#   Assembly and annotation of transcripts
############################################
# FLAIR Pipeline
hg38 = "/path/to/hg38_mainChr.fa"
genecodeV42 = "/path/to/gencode.v42.annotation.gtf"
flair = '/path/to/flair'
bam2bed = "/path/to/flair/bam2Bed12.py"
sort_bam = '/path/to/SampleID_ONT.sorted.bam'
sort_bed = '/path/to/flair_sorted.bed'
flair_out = '/path/to/assemble'
fastq = '/path/to/SampleID.ONT.fq'
correct_bed = '/path/to/flair_correct.bed'
temp_dir = '/path/to/flair_temp_dir'
sampleid_flair = '/path/to/flair_out/SampleID'
gffcompare = "/path/to/gffcompare"
flair_gtf = '/path/to/sampleid_flair/SampleID.isoform.gtf'

########### 1.align_bam2bed
'''
If your reads are already aligned, you can convert the sorted bam output to bed12 using bam2Bed12 to supply for flair-correct. This step smoothes gaps in the alignment.
'''
def flair_bam2bed(sort_bam, sort_bed):
    cmd = 'python {bam2bed} -i {sort_bam} > {sort_bed}'.format(bam2bed=bam2bed, sort_bam=sort_bam, sort_bed=sort_bed)
    os.system(cmd)

############# 2.correct
'''
This module corrects misaligned splice sites using genome annotations and/or short-read splice junctions.
'''
def flair_correct(sort_bed,flair_out):
    cmd = "{flair} correct -q {sort_bed} -g {hg38} -f {genecodeV42} -t 100 -o {flair_out}".format(
        flair=flair,sort_bed=sort_bed,hg38=hg38,genecodeV42=genecodeV42,flair_out=flair_out)
    os.system(cmd)

################ 3.collapse
def flair_collapse(flair,fastq,correct_bed,temp_dir,sampleid_flair):
    os.system('{flair} collapse -r {fastq} -q {correct_bed} -g {hg38} -f {gencode_gtf} -t 100 --generate_map '
              '--support 10 --quality 10 --stringent --temp_dir {temp_dir} -o {sampleid_flair}'.format(
        flair=flair, fastq=fastq, correct_bed=correct_bed, hg38=hg38,
        gencode_gtf=gtf, temp_dir=temp_dir, sampleid_flair=sampleid_flair))

################ 4.gffcompare
def transcript_annotated(output,flair_gtf):
    cmd = '{gffcompare} -R -r {genecodeV42} -o {output} {flair_gtf}'.format(gffcompare=gffcompare, genecodeV42=genecodeV42,output=output, flair_gtf=flair_gtf)
    os.system(cmd)


############################################
#            Identify TCTs
############################################
# coverage > 5
aligndir = "/path/to/align"
tedir = "/path/to/TE_read"
assembledir = "/path/to/assemble"
os.mkdir(tedir)

# get TE-related novel transcript bam
def novel_reads(sample, regionBed):
    os.chdir(tedir)
    samfile = f'/path/to/{sample}_ONT.sorted.bam'
    outbam = f'/path/to/{sample}_TE-related.novel-transcript.readsName.bam'
    os.system('samtools view -hb -L %s %s > %s' % (regionBed, samfile, outbam))
    os.system('samtools index %s_TE-related.novel-transcript.readsName.bam' % sample)


for dir in os.listdir(aligndir):
    print(dir)
    subdir = os.path.join(aligndir, dir)
    os.chdir(subdir)
    # TE-exon Reads
    te_reads = os.popen("awk -F '\t' '{if ($5>0) print $4}' %s | sort -u" % os.path.join(subdir,
                                                                                         'sorted_%s_ONT_rna-hg38TE.bed' % dir)).read().strip().split(
        "\n")
    exon_reads = os.popen("awk -F '\t' '{if ($5>0) print $4}' %s | sort -u" % os.path.join(subdir,
                                                                                           'sorted_%s_ONT_rna-exon.bed' % dir)).read().strip().split(
        "\n")
    common_reads = list(set(te_reads) & set(exon_reads))
    teReads=pd.read_csv('sorted_%s_ONT_rna-hg38TE.bed' % dir,sep='\t',header=None,usecols=[3,4,6,7,8,9,10,11,12])
    teReads=teReads[teReads[4]>10]
    teReads = teReads[teReads[3].isin(common_reads)]
    te_info = teReads[[6,7,8,9,10,11,12]]
    te_info['te'] = teReads[6].map(str)+'|'+teReads[7].map(str)+'|'+teReads[8].map(str)+'|'+teReads[10].map(str)+'|'+teReads[9].map(str)
    te_counts=te_info['te'].value_counts().rename_axis('unique_values').to_frame('counts').reset_index()
    te_counts = te_counts[te_counts['counts']>5]
    te_counts1 = te_counts['unique_values'].str.split('|',expand=True)
    te_counts1.to_csv(os.path.join(tedir,dir+'_ONT_rnaRead5-exon.hg38TE.bed'),sep='\t',index=False,header=False)
    # novel transcript
    tmap = pd.read_csv(os.path.join(assembledir,dir,'{sample}.{sample}.isoforms.gtf.tmap'.format(sample=dir)),sep='\t',header=0)
    tmap = tmap[tmap['num_exons']>1]
    tmap = tmap[~tmap['class_code'].isin(['=','x','i','e'])]
    gtf = os.path.join(assembledir,dir,'{sample}.isoforms.gtf'.format(sample=dir))
    newgtf = os.path.join(tedir,'{sample}.isoforms-novel.gtf'.format(sample=dir))
    qry_id = ' -e '.join(tmap['qry_id'])
    os.system("grep -e %s %s > %s" % (qry_id, gtf,newgtf))
    novel = pd.read_csv(newgtf,sep='\t',header=None)
    novel_exon = novel[novel[2]=='exon']
    novel_exon = novel_exon[[0,3,4,5,6,2,8]]
    novel_exon.to_csv(newgtf.split(".gtf")[0]+'-exon.bed',sep='\t',index=False,header=False,quoting=csv.QUOTE_NONE)
    os.system('bedtools intersect -a %s -b %s -wa -wb > %s' % (newgtf.split(".gtf")[0]+'-exon.bed',
                                                               os.path.join(tedir,dir+'_ONT_rnaRead5-exon.hg38TE.bed'),
                                                               newgtf.split(".gtf")[0]+'-exon.ONT_rnaRead5-exon.hg38TE.bed'))
    # TE-related novel transcript region bed
    tid = os.popen("awk -F 'transcript_id \"' '{print $2}' %s| awk -F '\"' '{print $1}'|sort -u " %
                   (newgtf.split(".gtf")[0]+'-exon.ONT_rnaRead5-exon.hg38TE.bed')).read().strip().split('\n')
    novel_transcript = novel[novel[2] == 'transcript']
    regionBed = os.path.join(tedir,dir+'_novelTranscript_ONT_rnaRead5-exon.hg38TE.bed')
    with open(regionBed,'w') as out:
        for i in range(0,novel_transcript.shape[0]):
            tid1 = novel_transcript.iloc[i,8].split('"')[3]
            if tid1 in tid:
                newline='\t'.join([novel_transcript.iloc[i,0],str(novel_transcript.iloc[i,3]-5000),
                                   str(novel_transcript.iloc[i,4]+5000),tid1,novel_transcript.iloc[i,8].split('"')[1]])+'\n'
                out.write(newline)
    novel_reads(dir, regionBed)
    sub_tedir = os.path.join(tedir, dir)
    if not os.path.exists(sub_tedir):
        os.mkdir(sub_tedir)
    os.chdir(sub_tedir)
    file_list_for_igv=['../%s_TE-related.novel-transcript.readsName.bam' % dir,
                       '../%s_TE-related.novel-transcript.readsName.bam.bai' % dir,
                       '../%s_ONT_rnaRead5-exon.hg38TE.bed' % dir,
                       newgtf,
                       newgtf.split(".gtf")[0]+'-exon.ONT_rnaRead5-exon.hg38TE.bed']
    for file in file_list_for_igv:
        os.system('mv %s ./' %file)

##################################################################
#      Generating a reference transcriptome including TCTs
##################################################################
os.mkdir(novel_assembledir)
novel = pd.read_csv(os.path.join(tedir,'5ONTreads-TE-novel-transcript.stat.txt'),header=0,sep='\t')
te_tid = os.popen("awk -F 'transcript_id \"' '{print $2}' %s/allCell.isoforms-novel-exon.ONT_rnaRead5-exon.hg38TE.bed \
                    | awk -F '\";' '{print $1}' | sort -u" % tedir).read().strip().split('\n')
novel = novel[novel['tid'].isin(te_tid)]

# get TE overlapped novel transcripts
for sample in os.listdir(aligndir):
    sub_novel = novel[novel['sample']==sample]
    gtf = os.path.join(assembledir,sample,sample+'.isoforms.gtf')
    new_gtf = os.path.join(novel_assembledir,sample+'.isoforms.5ONTreads-novelTE.gtf')
    tids = sub_novel['tid'].drop_duplicates()
    gtf_df = pd.read_csv(gtf, sep='\t',header=None)
    tid_list = pd.DataFrame(columns = gtf_df.columns)
    for tid in tids:
        tid_list=tid_list.append(gtf_df.loc[gtf_df[8].str.contains(tid)])
    tid_list.to_csv(new_gtf, sep='\t',header=False, index=False,quoting=csv.QUOTE_NONE)

os.system('stringtie --merge  \
            -G ../assemble/gencode.v42.annotation.gtf \
            -o gencodeV42.merged.allNovelTE.gtf \
            assembly_GTF_list.txt') # DRR228517.isoforms.5ONTreads-novelTE.gtf空

os.system('stringtie --merge  \
            -o merged.allNovelTE.gtf \
            new_assembly_GTF_list.txt') # DRR228517.isoforms.5ONTreads-novelTE.gtf空

merged_gtf=os.path.join(novel_assembledir,'merged.allNovelTE.gencode.v42.gtf')
new_merged_gtf = open(os.path.join(novel_assembledir,'addGeneName_merged.allNovelTE.gencode.v42.gtf'),'w')
Merge = open(merged_gtf)
line = Merge.readline()
gene_list = {}
while line:
    if 'MSTRG' in line:
        info = line.strip().split('\t')
        if info[2]=='transcript':
            gene_id = os.popen("grep '\<%s\>' %s/*.isoforms.5ONTreads-novelTE.gtf | grep -e '\<%s' -e '\<%s' \
                                | grep -v 'exon'" % (
                                info[0],novel_assembledir,info[3][0:6],info[4][0:6])
                                ).read().strip().split('\n')
            gene_ids = []
            for gene in gene_id:
                gene_ids.append(gene.split('gene_id "')[1].split('"')[0])
            gene_ids = list(set(gene_ids))
            if len(gene_ids)>1:
                print(info[0],info[3],info[4])
                break
            else:
                gene_id=gene_ids[0]
            if 'ENSG' in gene_ids[0]:
                gene_name = os.popen("grep '%s' %s/gencode.v42.annotation.gtf | head -1" %
                                    (gene_ids[0],assembledir)).read().strip().split('gene_name "')[1].split('"')[0]
            else:
                gene_name = gene_id
            oldid = info[-1].split('"')[1]
            gene_id = gene_ids[0]
            gene_list[oldid] = [gene_ids[0],gene_name]
        else:
            oldid = info[-1].split('"')[1]
            gene_id = gene_list[oldid][0]
            gene_name = gene_list[oldid][1]
        newline = line.split("; gene_id")[0]+'; gene_id"'+gene_id+'"; gene_name "'+gene_name+'"\n'
        new_merged_gtf.write(newline)
    else:
        new_merged_gtf.write(line)
    line = Merge.readline()



new_merged_gtf.close()
Merge.close()

##################################################################
#     Transcript-level quantification and candidate selection
##################################################################
# quantify
def quatification(sample):
    print(sample + ' start')
    subdir=os.path.join(novel_assembledir,sample)
    if not os.path.exists(subdir):
        os.mkdir(subdir)
    else:
        os.system('rm -r '+ subdir)
        os.mkdir(subdir)
    os.chdir(subdir)
    bam = os.path.join(aligndir,sample,'sorted_%s_ONT_rna.bam' % sample)
    os.system('stringtie -p 40 -L \
                -G {} \
                -o {}_merged.allNovelTE.gencode.v42.stringtie.gtf  \
                -b ballgown_out_dir -e \
                {}'.format(os.path.join(novel_assembledir,'addGeneName_merged.allNovelTE.gencode.v42.gtf'),
                           sample,sample, bam))
    print(sample+' end')


pools = Pool(3)
for sample in os.listdir(aligndir):
    pools.apply_async(quatification, args=(sample,))

pools.close()
pools.join()
del pools


# FPKM
for sample in os.listdir(aligndir):
    subdir=os.path.join(novel_assembledir,sample)
    os.chdir(subdir)
    tpm = pd.read_csv(os.path.join('ballgown_out_dir','t_data.ctab'),sep='\t',header=0)
    tpm1 = tpm[tpm['t_name'].str.contains('MSTRG')]
    novel_gene = list(set(tpm1['gene_name']))
    tpm2 = tpm[tpm['gene_name'].isin(novel_gene)]
    fraction_tpm = pd.DataFrame(columns=['sample','ENSG','genename','ENST','FPKM','fraction'])
    for gene in list(set(tpm2['gene_name'])):
        sub = tpm2[tpm2['gene_name']==gene]
        sub=sub[sub['FPKM']!=0]
        tids = sub['t_name']
        for tid in tids:
            f = sub[sub['t_name']==tid].iloc[0]['FPKM']/np.sum(sub['FPKM'])
            fraction_tpm.loc[len(fraction_tpm.index)] = [sample, gene,sub.iloc[0]['gene_id'],tid,
                                                         sub[sub['t_name']==tid].iloc[0]['FPKM'],str(f)]
    fraction_tpm.to_csv(sample+'_merged.allNovelTE.gencode.v42.stringtie.FPKM-Fraction.txt',sep='\t',header=False,index=False)
    print(sample+' end')

##################################################################
#         Open-reading frame prediction and annotation
##################################################################
TransDecoder_dir = "/path/to/transdecoder/bin"
LongOrfs = os.path.join(TransDecoder_dir,'TransDecoder.LongOrfs')
Predict = os.path.join(TransDecoder_dir,'TransDecoder.Predict')
gff3pl = os.path.join(TransDecoder_dir,'gtf_to_alignment_gff3.pl')
fastapl = os.path.join(TransDecoder_dir,'gtf_genome_to_cdna_fasta.pl')
ref_orf = os.path.join(TransDecoder_dir,'cdna_alignment_orf_to_genome_orf.pl')
gff3 = "/path/to/merged.allNovelTE.gff3"
transdecoder_gff3 = "/path/to/merged.allNovelTE.out.fa.transdecoder.gff3"
ref_gff3 = "/path/to/merged.allNovelTE.out.fa.transdecoder.genome.gff3"


def to_fa(gtf,hg38,outfa):
    cmd = '{fastapl} {gtf} {hg38} > {outfa}'.format(fastapl=fastapl,gtf=gtf,hg38=hg38,outfa=outfa)
    os.system(cmd)


def gff(gtf,gff3):
    cmd = '{gff3pl} {gtf} > {gff3}'.format(gff3pl=gff3pl,gtf=gtf,gff3=gff3)
    os.system(cmd)


def longOrf(LongOrfs,outfa):
    cmd = '{LongOrfs} -t {outfa}'.format(LongOrfs=LongOrfs,outfa=outfa)
    os.system(cmd)


def predict(Predict,outfa):
    cmd = '{Predict} -t {outfa}'.format(Predict=Predict,outfa=outfa)
    os.system(cmd)


def align(transdecoder_gff3,gff3,outfa,ref_gff3):
    cmd = '{ref_orf} {transdecoder_gff3} {gff3} {outfa} > {ref_gff3}'.format(ref_orf=ref_orf,transdecoder_gff3=transdecoder_gff3,
                                                                             gff3=gff3,outfa=outfa,ref_gff3=ref_gff3)
    os.system(cmd)
