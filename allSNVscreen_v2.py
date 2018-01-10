import os
import re
import sys
from time import time

import numpy as np
import pandas as pd
import xlsxwriter
from collections import defaultdict, Counter
sys.path.append('/home/wuj/bin/tools')
from processbar import printProgress


class GenoFromAllsnv:
    """Extracting `sample.geno` from allSNV.txt, and snvs are
    filtered by following rules tacitly:
        Chr          :   rm X, Y
        HomHit       :   1
        Indel        :   remove
        Region       :   re.match('exonic')
        CaseCallRate :   > 0.98
        CtlCallRate  :   > 0.98
        MAF          :   > 0.01
        HWE_ALL      :   > 0.0001
        HWE_Case     :   > 0.0001
        HWE_Ctl      :   > 0.0001
    """
    def __init__(self, allsnv,
                 sampleinfo,
                 Chrs=range(1,23),
                 homhit=1,
                 region='exonic',
                 case_callrate=0.98,
                 ctl_callrate=0.98,
                 maf=0.01,
                 hwe_all=0.0001,
                 hwe_case=0.0001,
                 hwe_ctl=0.0001,
                 outpath=None,
                 ):
        self.allsnv = allsnv
        self.sampleinfo = sampleinfo
        self.Chrs = Chrs
        self.homhit = homhit
        self.region = region
        self.case_callrate = case_callrate
        self.ctl_callrate = ctl_callrate
        self.maf = maf
        self.hwe_all = hwe_all
        self.hwe_case = hwe_case
        self.hwe_ctl = hwe_ctl
        self.outpath = outpath

        if outpath is None:
            self.outpath = './'

    def read_sampleinfo(self):
        dic = defaultdict(list)
        count = 0
        with open(self.sampleinfo, 'rt') as fh:
            for line in fh:
                count += 1
                if count == 1:
                    continue
                if re.match(r'^\s$', line):
                    continue
                arr = line.strip().split('\t')
                dic[arr[1]].append(arr[0])
        return dic

    def load_table_to_malformed(self):
        # table = pd.read_table(self.allsnv, header=0, index_col=0, sep='\t', nrows=10)
        allsnv_xlsx = os.path.join(self.outpath, 'allSNV-analysis.xlsx')
        workbook = xlsxwriter.Workbook(allsnv_xlsx)
        worksheet = workbook.add_worksheet()
        row_xlsx = 0

        anno_head = ['SNP ID', 'Chrs', 'Position', 'Ref Allele', 'Alt Allele', 'Gene']
        dic_sites = {}
        dic_genos = {}

        dic_samples = self.read_sampleinfo()

        total = os.path.getsize(self.allsnv)
        count = 0
        current_bytes_read = 0
        with open(self.allsnv, 'rt', encoding='utf-8', errors='ignore') as fh:
            for line in fh:
                current_bytes_read += len(line)
                printProgress(current_bytes_read, total, prefix='Progress: ', suffix='Complete', barLength=50, epsilon=0.000000001)
                count += 1
                line = re.sub(r'[\r\n]$','', line)
                arr = line.split('\t')
                if count == 1:
                    header = arr
                    samples_pos = header.index('Control Percentage') + 1
                    wanted_samples = dic_samples['case'] + dic_samples['control']
                    samples_index = [header.index(s) for s in wanted_samples]
                    wanted_cols_num = list(range(0, samples_pos)) + samples_index
                    wanted_cols_name = header[:samples_pos] + wanted_samples
                    row_xlsx = self.output_xlsx(wanted_cols_name, worksheet, row_xlsx)
                    # col_names_for_genodf = ['SNV NO.', 'Gene', 'Chrs', 'Position', 'Ref Allele', 'Alt Allele', 'SNP ID'] + wanted_samples
                    # geno_df = pd.DataFrame(columns=col_names_for_genodf)
                    continue
                line_series = pd.Series(self.item_by_index(arr, wanted_cols_num), index=wanted_cols_name)
                parsed_series = self.parse_line(line_series, dic_samples)
                if parsed_series is not None:
                    row_xlsx = self.output_xlsx(parsed_series, worksheet, row_xlsx)
                    # series_for_genodf = [parsed_series[name] for name in col_names_for_genodf]
                    # geno_df.loc[series_for_genodf[0]] = series_for_genodf
                    snp_id = parsed_series['SNP ID']
                    sites_info = '\t'.join([parsed_series[name] for name in anno_head])
                    dic_sites[snp_id] = sites_info
                    for sample in wanted_samples:
                        dic_genos.setdefault(sample, {})[snp_id] = parsed_series[sample].strip()
        self.output_geno(dic_sites, dic_genos, anno_head)
        workbook.close()

    @staticmethod
    def item_by_index(alist, index):
        return [alist[i] or '' for i in index]

    def output_xlsx(self, line, worksheet, row):
        for i, j in enumerate(line):
            worksheet.write(row, i, j)
        row += 1
        return row

    def output_geno(self, dic_sites, dic_genos, anno_head):
        fanno = os.path.join(self.outpath, 'anno.txt')
        fgeno = os.path.join(self.outpath, 'sample.geno')

        with open(fanno, 'wt') as fh:
            fh.write('\t'.join(anno_head) + '\n')
            for site in dic_sites:
                fh.write(dic_sites[site] + '\n')

        with open(fgeno, 'wt') as fh:
            all_sites = dic_sites.keys()
            fh.write('Samples\t' + '\t'.join(all_sites) + '\n')
            for sample in dic_genos:
                line = '\t'.join([dic_genos[sample].get(site) for site in all_sites])
                fh.write(sample + '\t' + line + '\n')

    def parse_line(self, series, dic_samples):
        genofmt = '{0}/{1}'
        cases = dic_samples['case']
        ctls = dic_samples['control']
        ref = series['Ref Allele']
        alt = series['Alt Allele']
        homr = genofmt.format(ref, ref)
        het = genofmt.format(ref, alt)
        homa = genofmt.format(alt, alt)
        case_geno = list(series[cases])
        ctl_geno = list(series[ctls])
        all_geno = case_geno + ctl_geno

        case_count = geno_count(case_geno, homr, het, homa)
        ctl_count = geno_count(ctl_geno, homr, het, homa)
        all_count = geno_count(all_geno, homr, het, homa)

        case_crate = callrate(case_geno)
        ctl_crate = callrate(ctl_geno)

        maf_all = maf_test(all_geno)
        hwe_all = HWE_exact_test(all_count[1], all_count[0], all_count[2])
        hwe_case = HWE_exact_test(case_count[1], case_count[0], case_count[2])
        hwe_ctl = HWE_exact_test(ctl_count[1], ctl_count[0], ctl_count[2])
        chrs = list(map(lambda x: str(x), self.Chrs))
        #if (series['Chrs'] in chrs and series['HOMOLOGY HITS'] == str(self.homhit) and
        #        re.search(self.region, series['Gene Region']) and not_indel(ref, alt) and
        #        maf_all > self.maf and case_crate > self.case_callrate and
        #        ctl_crate > self.ctl_callrate and hwe_all > self.hwe_all and
        #        hwe_case > self.hwe_case and hwe_ctl > self.hwe_ctl):
        if maf_all > 0.01:
        #if (maf_all > 0.05 and case_crate > 0.9 and ctl_crate > 0.9 and hwe_ctl > 0.001):

            rs_number = parse_rs(series['SNP ID'])
            if rs_number is None:
                series['SNP ID'] = series['SNV NO.']
            else:
                series['SNP ID'] = rs_number

            series['Case(0|1|2)'] = '|'.join(map(lambda x:str(x), case_count))
            series['Control(0|1|2)'] = '|'.join(map(lambda x:str(x), ctl_count))
            series['Genotyping Quality'] = 1
            alt_freq_all = alt_freq(all_geno, ref, alt)
            alt_freq_ctl = alt_freq(ctl_geno, ref, alt)
            series['Alt Allele Freq'] = alt_freq_all
            series['Control Percentage'] = alt_freq_ctl
            het_case, hom_case, het_ctl, hom_ctl = hethom_mapping(series, cases, ctls, het, homa)
            series['Case_het'] = het_case
            series['Case_hom'] = hom_case
            series['Control_het'] = het_ctl
            series['Control_hom'] = hom_ctl
            return series
        else:
            return None

def hethom_mapping(series, cases, ctls, het, homa):
    het_case = []
    het_ctl = []
    hom_case = []
    hom_ctl = []
    for s in cases:
        if series[s] == het:
            het_case.append(s)
        if series[s] == homa:
            hom_case.append(s)
    for s in ctls:
        if series[s] == het:
            het_ctl.append(s)
        if series[s] == homa:
            hom_ctl.append(s)
    return ','.join(het_case), ','.join(hom_case), ','.join(het_ctl), ','.join(hom_ctl)

def not_indel(ref, alt):
    tmp_str = ref + alt
    if re.search(r'-', tmp_str):
        return False
    return True

def parse_rs(snvid):
    rs = re.match(r'rs\d+', snvid)
    try:
        return rs.group()
    except:
        return None

def geno_count(genos, homr, het, homa):
    homr_num = genos.count(homr)
    het_num = genos.count(het)
    homa_num = genos.count(homa)
    return [homr_num, het_num, homa_num]

def callrate(genos):
    rest = list(filter(lambda x: x.strip(), genos))
    try:
        rate = len(rest)/len(genos)
    except ZeroDivisionError:
        rate = 0
    return rate

def maf_test(genos):
    try:
        alleles = []
        for geno in genos:
            arr = geno.split('/')
            for a in arr:
                if a.strip():
                    alleles.append(a)
        count = Counter(alleles)
        count_sorted = sorted(count.items(), key=lambda l: l[1])
        if len(count) < 2:
            return 0
        return count_sorted[0][1]/len(alleles)
    except:
        return 0

def alt_freq(genos, ref, alt):
    try:
        sum_all = 0
        sum_alt = 0
        for geno in genos:
            arr = geno.split('/')
            for a in arr:
                if a.strip():
                    sum_all += 1
                if a == alt:
                    sum_alt += 1
        return sum_alt/sum_all
    except:
        return 0

def HWE_exact_test(obs_het, obs_hom1, obs_hom2):
    try:
        obs_homr = obs_hom2 if obs_hom1 < obs_hom2 else obs_hom1
        obs_homa = obs_hom1 if obs_hom1 < obs_hom2 else obs_hom2

        rare_copies = 2 * obs_homa + obs_het
        all_genotypes = obs_het + obs_homa + obs_homr

        het_probs = [0.0] * (rare_copies + 1)

        # start at midpoint
        mid = int(rare_copies * (2 * all_genotypes - rare_copies)/(2 * all_genotypes))

        # check to ensure that midpoint and rare alleles have same parity
        if (rare_copies & 1) ^ (mid & 1):
            mid += 1

        curr_het = mid
        curr_homa = (rare_copies - mid)/2
        curr_homr = all_genotypes - curr_het - curr_homa

        het_probs[mid] = 1.0
        sum_ = float(het_probs[mid])

        for curr_het in range(mid, 1, -2):
            het_probs[curr_het - 2] = het_probs[curr_het] * curr_het * (curr_het - 1.0)/(4 * (curr_homa + 1.0) * (curr_homr + 1.0))
            sum_ += het_probs[curr_het - 2]

            # 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
            curr_homa += 1
            curr_homr += 1

        curr_het = mid
        curr_homa = (rare_copies - mid)/2
        curr_homr = all_genotypes - curr_het - curr_homa

        for curr_het in range(mid, rare_copies - 1, 2):
            het_probs[curr_het + 2] = het_probs[curr_het] * 4.0 * curr_homa * curr_homr/((curr_het + 2.0) * (curr_het + 1.0))
            sum_ += het_probs[curr_het + 2]

            # add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
            curr_homr -= 1
            curr_homa -= 1

        for i in range(0, rare_copies + 1):
            het_probs[i] /= sum_

        # alternate p-value calculation for p_hi/p_lo
        p_hi = float(het_probs[obs_het])
        for i in range(obs_het, rare_copies + 1):
            p_hi += het_probs[i]

        p_lo = float(het_probs[obs_het])
        for i in range(obs_het - 1, -1, -1):
            p_lo += het_probs[i]

        p_hi_lo = 2.0 * p_hi if p_hi < p_lo else 2.0 * p_lo

        p_hwe = 0.0
        # p-value calculation for p_hwe
        for i in range(0, rare_copies + 1):
            if het_probs[i] > het_probs[obs_het]:
                continue
            p_hwe += het_probs[i]

        p_hwe = 1.0 if p_hwe > 1.0 else p_hwe
        return p_hwe
    except:
        return 0



if __name__ == '__main__':
    t1 = time()
    sampleinfo, allsnv = sys.argv[1:]
    screener = GenoFromAllsnv(allsnv, sampleinfo)
    screener.load_table_to_malformed()
    print('time cost:', time() - t1)


