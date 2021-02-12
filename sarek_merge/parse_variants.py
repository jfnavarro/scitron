"""
A script that parses annotaded variants from Sarek and generates an unique matrix
with the variants found per sample. The table also contains extra columns with
useful information.
"""

from collections import defaultdict
from pathlib import Path
import pandas as pd
from variants import *
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import subprocess
import sys
from itertools import chain


def main(FILTER_DP,
         FILTER_DP_PDX,
         FILTER_VAF,
         SAMPLES,
         PATH,
         DISABLE_EFFECT_FILTER):

    metasamples = pd.read_csv(SAMPLES, sep="\t", header=None)
    metasamples.columns = ['PATIENT_ID', 'GENDER', "TUMOR", "SAMPLE_ID", "LANE", "FASTQ1", "FASTQ2"]
    metasamples = metasamples.drop_duplicates(subset=['SAMPLE_ID'])

    print('Loaded {} samples'.format(metasamples.shape[0]))
    print('Searching for variants in {}'.format(PATH))

    # FIRST parse VCF files and create a dictionary by variant key containing the variants info
    variants_dict = defaultdict(list)
    for index, row in metasamples.iterrows():
        sample_id = row['SAMPLE_ID']
        sample_patient = row['PATIENT_ID']
        tumor_only = all(metasamples[metasamples['PATIENT_ID'] == sample_patient]['TUMOR'] == 1)
        germline = all(metasamples[metasamples['PATIENT_ID'] == sample_patient]['TUMOR'] == 0)
        # Retrieve variants
        if tumor_only or germline:
            # TODO make this a single reg exp
            for elem in chain(Path(PATH).rglob('*{}_snpEff.ann.vcf.gz'.format(sample_id)),
                              Path(PATH).rglob('*{}_variants_snpEff.ann.vcf.gz'.format(sample_id))):
                variants_list = list()
                if 'Strelka' in str(elem):
                    print('Parsing variant in {}'.format(str(elem)))
                    variants_list = parse_strelka_germline_variants(elem, FILTER_DP, FILTER_DP_PDX,
                                                                    FILTER_VAF, DISABLE_EFFECT_FILTER)
                if 'HaplotypeCaller' in str(elem) and '.g_snpEff' not in str(elem):
                    print('Parsing variant in {}'.format(str(elem)))
                    variants_list = parse_haplotypecaller_germline_variants(elem, FILTER_DP, FILTER_DP_PDX,
                                                                            FILTER_VAF, DISABLE_EFFECT_FILTER)
                if len(variants_list) > 0:
                    for v in variants_list:
                        v.sample = sample_id
                        v.patient = sample_patient
                        variants_dict[v.key].append(v)
        else:
            gdna = metasamples[(metasamples['PATIENT_ID'] == sample_patient)
                               & (metasamples['TUMOR'] == 0)]['SAMPLE_ID'].to_numpy()[0]
            #gdna = gdna[:gdna.rfind('_')]
            if row['TUMOR'] != 0:
                # TODO make this a single reg exp
                for elem in chain(Path(PATH).rglob('*{}_vs_{}_snpEff.ann.vcf.gz'.format(sample_id, gdna)),
                                  Path(PATH).rglob('*{}_vs_{}_variants_snpEff.ann.vcf.gz'.format(sample_id, gdna))):
                    variants_list = list()
                    if 'StrelkaBP' in str(elem) and 'somatic_snvs' in str(elem):
                        print('Parsing variant in {}'.format(str(elem)))
                        variants_list = parse_strelka_variants(elem, FILTER_DP, FILTER_DP_PDX,
                                                               FILTER_VAF, DISABLE_EFFECT_FILTER)
                    if 'StrelkaBP' in str(elem) and 'somatic_indels' in str(elem):
                        print('Parsing variant in {}'.format(str(elem)))
                        variants_list = parse_strelka_indel_variants(elem, FILTER_DP, FILTER_DP_PDX,
                                                                     FILTER_VAF, DISABLE_EFFECT_FILTER)
                    if 'Mutect2_filtered' in str(elem):
                        print('Parsing variant in {}'.format(str(elem)))
                        variants_list = parse_mutect_variants(elem, FILTER_DP, FILTER_DP_PDX,
                                                              FILTER_VAF, DISABLE_EFFECT_FILTER)
                    if len(variants_list) > 0:
                        for v in variants_list:
                            v.sample = sample_id
                            v.patient = sample_patient
                            variants_dict[v.key].append(v)

    # SECOND parse the dictionary to create a table with the variants and useful information
    ALL_SAMPLES = list(metasamples['SAMPLE_ID'])
    ALL_PATIENTS = list(metasamples['PATIENT_ID'])
    with open('merged_variants.txt', 'w') as handler:
        with open('tmp_anno.txt', 'w') as handler2:
            header_samples = '\t'.join(['{}-{}'.format(x, y) for x, y in zip(ALL_SAMPLES, ALL_PATIENTS)])
            handler.write(
                'VARIANT_KEY\tGENE\tEFFECT\tAACHANGE\tTYPE\tNUM_SAMPLES\tNUM_PATIENTS\t{}\n'.format(header_samples))
            for key, items in variants_dict.items():
                samples_variants_filtered = defaultdict(list)
                effects = set()
                patients = set()
                for v in items:
                    effects.update(v.effects)
                    patients.add(v.patient)
                    samples_variants_filtered[v.sample].append('{}:{};{};{}'.format(v.caller, v.DP, v.AD, v.VAF))
                if len(samples_variants_filtered) >= 1:
                    key_split = key.split(':')
                    chr = key_split[0]
                    pos = key_split[1].split()[0]
                    ref = key_split[1].split()[1].split('>')[0]
                    alt = key_split[1].split()[1].split('>')[1]
                    num_patients = len(patients)
                    num_samples = len(samples_variants_filtered)
                    samples_str = '\t'.join(['|'.join(set(samples_variants_filtered[s]))
                                             if s in samples_variants_filtered else 'Na' for s in ALL_SAMPLES])
                    to_write = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(key,
                                                                         '|'.join(set([x[1] for x in effects])),
                                                                         '|'.join(set([x[0] for x in effects])),
                                                                         '|'.join(set([x[2] for x in effects])),
                                                                         '|'.join(set([x[3] for x in effects])),
                                                                         num_samples,
                                                                         num_patients,
                                                                         samples_str)
                    handler.write(to_write)
                    handler2.write('{}\t{}\t{}\t{}\t{}\n'.format(chr, pos, pos, ref, alt))

    # Annotate with the annovar the variants
    ANNOVAR_PATH = os.environ['ANNOVAR_PATH']
    INPUT = 'tmp_anno.txt'
    THREADS = 4

    cmd = '{} {} {} -buildver hg38 -thread {} -out ann ' \
          '-protocol avsnp150,gnomad_exome,cosmic70,clinvar_20200316,dbnsfp41c -operation f,f,f,f,f -nastring Na'.format(
        os.path.join(ANNOVAR_PATH, 'table_annovar.pl'),
        INPUT,
        os.path.join(ANNOVAR_PATH, 'humandb'),
        THREADS)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, error = p.communicate()
    if p.returncode != 0:
        print(output)
        print(error)
        sys.exit(-1)

    # Load Annovar variants info
    anno_dict = dict()
    with open('ann.hg38_multianno.txt', 'r') as handler:
        lines = handler.readlines()
        header = lines.pop(0).rstrip().split('\t')
        for line in lines:
            tokens = line.split()
            key = '{}:{} {}>{}'.format(tokens[header.index('Chr')],
                                       tokens[header.index('Start')],
                                       tokens[header.index('Ref')],
                                       tokens[header.index('Alt')])
            value = (tokens[header.index('avsnp150')],
                     tokens[header.index('gnomAD_exome_ALL')],
                     tokens[header.index('cosmic70')],
                     '{}:{}'.format(tokens[header.index('CLNSIG')],
                                    tokens[header.index('CLNDISDB')]),
                     tokens[header.index('MutationAssessor_pred')],
                     tokens[header.index('MutationTaster_pred')])
            anno_dict[key] = value

    # Add info to table
    variants_table = pd.read_csv('merged_variants.txt', sep="\t", index_col=0)
    gene_patients = defaultdict(list)
    gene_samples = defaultdict(list)
    for index, row in variants_table.iterrows():
        for c in variants_table.columns:
            if row[c] != 'Na':
                gene_patients[row['GENE']].append(c.split('-')[-1])
                gene_samples[row['GENE']].append(c)
    variants_table.insert(4, 'DBSNP', [anno_dict[v][0] for v in variants_table.index])
    variants_table.insert(5, 'GNOMAD', [anno_dict[v][1] for v in variants_table.index])
    variants_table.insert(6, 'COSMI70', [anno_dict[v][2] for v in variants_table.index])
    variants_table.insert(7, 'CLINVAR', [anno_dict[v][3] for v in variants_table.index])
    variants_table.insert(8, 'MutationAssessor_pred', [anno_dict[v][4] for v in variants_table.index])
    variants_table.insert(9, 'MutationTaster_pred', [anno_dict[v][5] for v in variants_table.index])
    variants_table.insert(11, 'NUM_SAMPLES_GENE', [len(set(gene_samples[v])) for v in variants_table['GENE']])
    variants_table.insert(13, 'NUM_PATIENTS_GENE', [len(set(gene_patients[v])) for v in variants_table['GENE']])

    variants_table.to_csv('merged_variants_anno.txt', sep="\t")


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--filter-dp', type=int, default=50, required=False,
                        help='Value for the DP filter (coverage). Default=50')
    parser.add_argument('--filter-dp-pdx', type=int, default=25, required=False,
                        help='Value for the DP filter for PDX samples if any (coverage). Default=25')
    parser.add_argument('--filter-vaf', type=float, default=5, required=False,
                        help='Value for the VAF filter (Variant Allele Frequency). Default=5')
    parser.add_argument('--samples', type=str, default=None, required=True,
                        help='Path to the file with the samples info (metadata Sarek).')
    parser.add_argument('--path', type=str, default=None, required=True,
                        help='Path to the folder of the results from Sarek (Annotated folder).')
    parser.add_argument('--disable-effect-filter', action='store_true', default=False, required=False,
                        help='Disable the effect filter that includes only missense variants and frameshifts.')
    args = parser.parse_args()
    main(args.filter_dp,
         args.filter_dp_pdx,
         args.filter_vaf,
         os.path.abspath(args.samples),
         os.path.abspath(args.path),
         args.disable_effect_filter)