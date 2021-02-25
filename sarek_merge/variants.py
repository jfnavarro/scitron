import os
import vcfpy
import numpy as np


class Variant:

    def __init__(self):
        self.chrom = None
        self.start = None
        self.ref = None
        self.alt = None
        self.effects = None
        self.caller = None
        self.VAF = None
        self.DP = None
        self.AD = None
        self.sample = None
        self.patient = None
        self.status = None

    @property
    def key(self):
        return "{}:{} {}>{}".format(self.chrom, self.start, self.ref, self.alt)

    def __str__(self):
        return '{}:{} {}>{}'.format(self.chrom, self.start, self.ref, self.alt)


def parse_haplotypecaller_germline_variants(filename, FILTER_DP, FILTER_DP_PDX, FILTER_VAF, NO_EFFECT_FILTER=False):
    variants = list()
    tumor_id = str(os.path.basename(filename)).split('_snpEff')[0].split('HaplotypeCaller_')[1]
    reader = vcfpy.Reader.from_path(filename)
    FDP = FILTER_DP_PDX if 'PDX' in str(filename) else FILTER_DP
    for record in reader:
        called = {x.sample: x.data for x in record.calls}
        tumor_DP = int(called[tumor_id]['DP'])
        tumor_AD = int(called[tumor_id]['AD'][1])
        tumor_VAF = np.around(tumor_AD / float(tumor_DP) * 100, 3) if tumor_DP > 0.0 else 0.0
        if tumor_VAF >= FILTER_VAF and tumor_DP >= FDP:
            effects = set()
            for effect in record.INFO['ANN']:
                gene = effect.split('|')[3]
                effect_name = effect.split('|')[1]
                effect_type = effect.split('|')[2]
                effect_biotype = effect.split('|')[7]
                effect_aachange = effect.split('|')[10]
                effect_feature = effect.split('|')[5]
                if NO_EFFECT_FILTER or effect_type in ['HIGH', 'MODERATE']:
                    effects.add((effect_name, gene, effect_aachange, effect_feature))
            if len(effects) > 0:
                variant = Variant()
                variant.chrom = record.CHROM
                variant.start = record.POS
                variant.ref = record.REF
                variant.alt = record.ALT[0].serialize()
                variant.caller = 'HaplotypeCaller'
                variant.effects = effects
                variant.DP = tumor_DP
                variant.AD = tumor_AD
                variant.VAF = tumor_VAF
                variants.append(variant)
    return variants


def parse_strelka_germline_variants(filename, FILTER_DP, FILTER_DP_PDX, FILTER_VAF, NO_EFFECT_FILTER=False):
    variants = list()
    tumor_id = str(os.path.basename(filename)).split('_variants_snpEff')[0].split('Strelka_')[1]
    reader = vcfpy.Reader.from_path(filename)
    FDP = FILTER_DP_PDX if 'PDX' in str(filename) else FILTER_DP
    for record in reader:
        called = {x.sample: x.data for x in record.calls}
        DP_FIELD = 'DP' if 'SNV' in record.ALT[0].type else 'DPI'
        tumor_DP = int(called[tumor_id][DP_FIELD])
        tumor_AD = int(called[tumor_id]['AD'][1])
        tumor_VAF = np.around((tumor_AD / float(tumor_DP)) * 100, 3) if tumor_DP > 0.0 else 0.0
        if 'PASS' in record.FILTER and tumor_VAF >= FILTER_VAF and tumor_DP >= FDP:
            effects = set()
            for effect in record.INFO['ANN']:
                gene = effect.split('|')[3]
                effect_name = effect.split('|')[1]
                effect_type = effect.split('|')[2]
                effect_biotype = effect.split('|')[7]
                effect_aachange = effect.split('|')[10]
                effect_feature = effect.split('|')[5]
                if NO_EFFECT_FILTER or effect_type in ['HIGH', 'MODERATE']:
                    effects.add((effect_name, gene, effect_aachange, effect_feature))
            if len(effects) > 0:
                variant = Variant()
                variant.chrom = record.CHROM
                variant.start = record.POS
                variant.ref = record.REF
                variant.alt = record.ALT[0].serialize()
                variant.caller = 'StrelkaGermline'
                variant.effects = effects
                variant.DP = tumor_DP
                variant.AD = tumor_AD
                variant.VAF = tumor_VAF
                variants.append(variant)
    return variants


def parse_strelka_variants(filename, FILTER_DP, FILTER_DP_PDX, FILTER_VAF, NO_EFFECT_FILTER=False):
    variants = list()
    reader = vcfpy.Reader.from_path(filename)
    FDP = FILTER_DP_PDX if 'PDX' in str(filename) else FILTER_DP
    for record in reader:
        called = {x.sample: x.data for x in record.calls}
        ref_index = record.REF + 'U'
        alt_index = str(record.ALT[0].serialize()) + 'U'
        normal_AD1 = int(called['NORMAL'][ref_index][0])
        normal_AD2 = int(called['NORMAL'][alt_index][0])
        normal_DP = normal_AD1 + normal_AD2
        normal_VAF = np.around((normal_AD2 / float(normal_DP)) * 100, 3) if normal_DP > 0.0 else 0.0
        tumor_AD1 = int(called['TUMOR'][ref_index][0])
        tumor_AD2 = int(called['TUMOR'][alt_index][0])
        tumor_DP = tumor_AD1 + tumor_AD2
        tumor_VAF = np.around((tumor_AD2 / float(tumor_DP)) * 100, 3) if tumor_DP > 0.0 else 0.0
        if 'PASS' in record.FILTER and 'SOMATIC' in record.INFO \
                and tumor_VAF >= FILTER_VAF and tumor_DP >= FDP:
            effects = set()
            for effect in record.INFO['ANN']:
                gene = effect.split('|')[3]
                effect_name = effect.split('|')[1]
                effect_type = effect.split('|')[2]
                effect_biotype = effect.split('|')[7]
                effect_aachange = effect.split('|')[10]
                effect_feature = effect.split('|')[5]
                if NO_EFFECT_FILTER or effect_type in ['HIGH', 'MODERATE']:
                    effects.add((effect_name, gene, effect_aachange, effect_feature))
            if len(effects) > 0:
                variant = Variant()
                variant.chrom = record.CHROM
                variant.start = record.POS
                variant.ref = record.REF
                variant.alt = record.ALT[0].serialize()
                variant.caller = 'Strelka'
                variant.effects = effects
                variant.DP = tumor_DP
                variant.AD = tumor_AD2
                variant.VAF = tumor_VAF
                variants.append(variant)
    return variants


def parse_strelka_indel_variants(filename, FILTER_DP, FILTER_DP_PDX, FILTER_VAF, NO_EFFECT_FILTER=False):
    variants = list()
    reader = vcfpy.Reader.from_path(filename)
    FDP = FILTER_DP_PDX if 'PDX' in str(filename) else FILTER_DP
    for record in reader:
        called = {x.sample: x.data for x in record.calls}
        normal_AD1 = int(called['NORMAL']['TAR'][0])
        normal_AD2 = int(called['NORMAL']['TIR'][0])
        normal_DP = normal_AD1 + normal_AD2
        normal_VAF = np.around((normal_AD2 / float(normal_DP)) * 100, 3) if normal_DP > 0.0 else 0.0
        tumor_AD1 = int(called['TUMOR']['TAR'][0])
        tumor_AD2 = int(called['TUMOR']['TIR'][0])
        tumor_DP = tumor_AD1 + tumor_AD2
        tumor_VAF = np.around((tumor_AD2 / float(tumor_DP)) * 100, 3) if tumor_DP > 0.0 else 0.0
        if 'PASS' in record.FILTER and 'SOMATIC' in record.INFO \
                and tumor_VAF >= FILTER_VAF and tumor_DP >= FDP:
            effects = set()
            for effect in record.INFO['ANN']:
                gene = effect.split('|')[3]
                effect_name = effect.split('|')[1]
                effect_type = effect.split('|')[2]
                effect_biotype = effect.split('|')[7]
                effect_aachange = effect.split('|')[10]
                effect_feature = effect.split('|')[5]
                if NO_EFFECT_FILTER or effect_type in ['HIGH', 'MODERATE']:
                    effects.add((effect_name, gene, effect_aachange, effect_feature))
            if len(effects) > 0:
                variant = Variant()
                variant.chrom = record.CHROM
                variant.start = record.POS
                variant.ref = record.REF
                variant.alt = record.ALT[0].serialize()
                variant.caller = 'Strelka-indel'
                variant.effects = effects
                variant.DP = tumor_DP
                variant.AD = tumor_AD2
                variant.VAF = tumor_VAF
                variants.append(variant)
    return variants


def parse_mutect_variants(filename, FILTER_DP, FILTER_DP_PDX, FILTER_VAF, NO_EFFECT_FILTER=False):
    variants = list()
    tumor_id = str(os.path.basename(filename)).split('_vs_')[0].split('Mutect2_filtered_')[1]
    normal_id = str(os.path.basename(filename)).split('_vs_')[1].split('_snpEff')[0]
    FDP = FILTER_DP_PDX if 'PDX' in str(filename) else FILTER_DP
    reader = vcfpy.Reader.from_path(filename)
    for record in reader:
        called = {x.sample: x.data for x in record.calls}
        normal_DP = int(called[normal_id]['DP'])
        normal_AD = int(called[normal_id]['AD'][1])
        normal_VAF = np.around(float(called[normal_id]['AF'][0]) * 100, 3) if normal_DP > 0.0 else 0.0
        tumor_DP = int(called[tumor_id]['DP'])
        tumor_AD = int(called[tumor_id]['AD'][1])
        tumor_VAF = np.around(float(called[tumor_id]['AF'][0]) * 100, 3) if tumor_DP > 0.0 else 0.0
        if 'PASS' in record.FILTER and tumor_VAF >= FILTER_VAF and tumor_DP >= FDP:
            effects = set()
            for effect in record.INFO['ANN']:
                gene = effect.split('|')[3]
                effect_name = effect.split('|')[1]
                effect_type = effect.split('|')[2]
                effect_biotype = effect.split('|')[7]
                effect_aachange = effect.split('|')[10]
                effect_feature = effect.split('|')[5]
                if NO_EFFECT_FILTER or effect_type in ['HIGH', 'MODERATE']:
                    effects.add((effect_name, gene, effect_aachange, effect_feature))
            if len(effects) > 0:
                variant = Variant()
                variant.chrom = record.CHROM
                variant.start = record.POS
                variant.ref = record.REF
                variant.alt = record.ALT[0].serialize()
                variant.caller = 'Mutect'
                variant.effects = effects
                variant.DP = tumor_DP
                variant.AD = tumor_AD
                variant.VAF = tumor_VAF
                variants.append(variant)
    return variants


def parse_mutect_tumor_onnly_variants(filename, FILTER_DP, FILTER_DP_PDX, FILTER_VAF, NO_EFFECT_FILTER=False):
    variants = list()
    tumor_id = str(os.path.basename(filename)).split('Mutect2_filtered_')[1].replace('_snpEff.ann.vcf.gz', '')
    FDP = FILTER_DP_PDX if 'PDX' in str(filename) else FILTER_DP
    reader = vcfpy.Reader.from_path(filename)
    for record in reader:
        called = {x.sample: x.data for x in record.calls}
        tumor_DP = int(called[tumor_id]['DP'])
        tumor_AD = int(called[tumor_id]['AD'][1])
        tumor_VAF = np.around(float(called[tumor_id]['AF'][0]) * 100, 3) if tumor_DP > 0.0 else 0.0
        if 'PASS' in record.FILTER and tumor_VAF >= FILTER_VAF and tumor_DP >= FDP:
            effects = set()
            for effect in record.INFO['ANN']:
                gene = effect.split('|')[3]
                effect_name = effect.split('|')[1]
                effect_type = effect.split('|')[2]
                effect_biotype = effect.split('|')[7]
                effect_aachange = effect.split('|')[10]
                effect_feature = effect.split('|')[5]
                if NO_EFFECT_FILTER or effect_type in ['HIGH', 'MODERATE']:
                    effects.add((effect_name, gene, effect_aachange, effect_feature))
            if len(effects) > 0:
                variant = Variant()
                variant.chrom = record.CHROM
                variant.start = record.POS
                variant.ref = record.REF
                variant.alt = record.ALT[0].serialize()
                variant.caller = 'MutectTumorOnly'
                variant.effects = effects
                variant.DP = tumor_DP
                variant.AD = tumor_AD
                variant.VAF = tumor_VAF
                variants.append(variant)
    return variants