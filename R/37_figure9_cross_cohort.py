#!/usr/bin/env python3
# 37_figure9_cross_cohort.py —— 跨队列验证合成图（Figure 8，去 trap 标签版）
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess, os

BASE = os.path.dirname(os.path.abspath(__file__))
F_21653 = os.path.join(BASE, 'GSE21653_series_matrix.txt')
F_ANN   = os.path.join(BASE, 'GPL570-55999.txt')
F_25066 = os.path.join(BASE, '35_gse25066_trap_assignment.csv')
OUT     = os.path.join(BASE, 'Figure_8_cross_cohort_validation.png')

# ---------- GSE21653 解析（与 36 号分析同口径：全队列 z 标准化） ----------
n0 = int(subprocess.run(['grep','-n','!series_matrix_table_begin',F_21653],
        capture_output=True, text=True).stdout.split(':')[0])
n1 = int(subprocess.run(['grep','-n','!series_matrix_table_end',F_21653],
        capture_output=True, text=True).stdout.split(':')[0])
dfm = pd.read_csv(F_21653, sep='\t', skiprows=n0, nrows=n1-n0-1, header=None,
                  low_memory=False, na_values=['','NA','"NA"'])
probe_ids = dfm.iloc[:,0].astype(str).values
mat = dfm.iloc[:,1:].apply(pd.to_numeric, errors='coerce').values
if probe_ids[0] in ('ID_REF','ID') or str(probe_ids[0]).startswith('"') or 'GSM' in str(probe_ids[1]):
    mat = mat[1:]; probe_ids = probe_ids[1:]
lines = [l[:-1] if l.endswith('\r') else l for l in open(F_21653, encoding='utf-8', errors='replace').read().split('\n')]
gsms = [x.strip('"') for x in [l for l in lines if l.startswith('!Sample_geo_accession')][0].split('\t')[1:]]
expr = pd.DataFrame(mat, index=probe_ids, columns=gsms)
fields = [dict() for _ in range(expr.shape[1])]
for cl in [l for l in lines if l.startswith('!Sample_characteristics_ch1')]:
    for j, v in enumerate([x.strip('"') for x in cl.split('\t')[1:]]):
        if ': ' in v:
            k, val = v.split(': ', 1); fields[j][k.strip()] = val.strip()
clin = pd.DataFrame(fields, index=gsms)
ann = pd.read_csv(F_ANN, sep='\t', comment='#', dtype=str)[['ID','Gene Symbol']].dropna()
ann = ann[ann['Gene Symbol'].str.strip() != '']
ann['symbol'] = ann['Gene Symbol'].str.split(' /// ').str[0].str.strip()
ann = ann[ann['ID'].isin(expr.index)]
ann2 = ann.copy(); ann2['em'] = ann2['ID'].map(expr.mean(axis=1))
idx_best = ann2.sort_values('em', ascending=False).drop_duplicates('symbol').set_index('symbol')['ID']
e = expr.loc[list(idx_best.values)]; e.index = list(idx_best.index)
tnbc = clin.index[(clin['er ihc']=='0') & (clin['pr ihc']=='0') & (clin['erbb2']=='0')]
ez = e.sub(e.mean(axis=1), axis=0).div(e.std(axis=1).replace(0,np.nan), axis=0)
tscore = ez.loc[['CD8A','CD8B','CD3D','CD3E','CD4','IL7R'], tnbc].mean(axis=0)
ido1z  = ez.loc['IDO1', tnbc]
d53 = pd.DataFrame({'t': pd.to_numeric(clin.loc[tnbc,'dfs time (months)'], errors='coerce')/12,
                    'e': pd.to_numeric(clin.loc[tnbc,'dfs evt'], errors='coerce'),
                    'trap': ((tscore<0)&(ido1z<0)).values}, index=tnbc).dropna()
assert len(d53)==85 and int(d53.trap.sum())==24, 'GSE21653 复现校验失败'

d66 = pd.read_csv(F_25066)
d66['sample'] = d66['sample'].astype(str).str.replace('"','', regex=False)
d66 = d66.dropna(subset=['drfs_time'])
assert len(d66)==117

# ---------- 绘图 ----------
def km(t, e_, ax, color, lab, lw=1.8):
    ts = np.concatenate([[0], np.sort(np.unique(t[e_==1]))]); s=1.0; S=[1.0]; xs=[0]
    for x in ts[1:]:
        n=(t>=x).sum(); dd=e_[t==x].sum()
        if n>0: s*=max(0,1-dd/n)
        S.append(s); xs.append(x)
    ax.step(xs, S, where='post', color=color, label=lab, lw=lw)

fig = plt.figure(figsize=(14, 10.5))
gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.26)

axA = fig.add_subplot(gs[0,0])
coh = ['TCGA C1\n(n=123)','METABRIC M1\n(n=320)','GSE21653\n(n=85)','I-SPY2\n(n=114 TNBC)','GSE25066\n(n=182)']
pct = [23, 30, 28.2, 25.4, 34.6]
axA.bar(range(5), pct, color=['#4C72B0','#4C72B0','#55A868','#DD8452','#C44E52'])
for i,v in enumerate(pct): axA.text(i, v+0.8, f'{v}%', ha='center', fontsize=9.5, fontweight='bold')
axA.axhspan(23, 35, color='grey', alpha=0.12)
axA.set_xticks(range(5)); axA.set_xticklabels(coh, fontsize=8.5)
axA.set_ylabel('CD8-low/IDO1-low profile (%)'); axA.set_ylim(0, 46)
axA.set_title('A  CD8-low/IDO1-low profile prevalence across five cohorts',
              loc='left', fontsize=11.5, fontweight='bold')

axB = fig.add_subplot(gs[0,1])
for g,c,l in [(True,'#E41A1C',f'CD8-low/IDO1-low (n={int(d53.trap.sum())})'),
              (False,'#377EB8',f'Others (n={int((~d53.trap).sum())})')]:
    m = d53.trap.values==g
    km(d53.t.values[m], d53.e.values[m].astype(int), axB, c, l)
axB.set_title('B  GSE21653: disease-free survival', loc='left', fontsize=11.5, fontweight='bold')
axB.set_xlabel('Time (years)'); axB.set_ylabel('DFS probability'); axB.set_ylim(0,1.03)
axB.legend(fontsize=9, loc='lower left')
axB.text(0.97, 0.62, 'log-rank p < 0.001\nHR 4.84 (95% CI 2.23-10.50)\nn = 85, 27 events',
         transform=axB.transAxes, ha='right', fontsize=9,
         bbox=dict(boxstyle='round', fc='white', ec='grey', alpha=0.9))

axC1 = fig.add_subplot(gs[1,0])
axC1.bar([0,1],[25.6,38.2], color=['#E41A1C','#377EB8'], width=0.55)
axC1.text(0,27.5,'25.6%\n(10/39)', ha='center', fontsize=9)
axC1.text(1,40.5,'38.2%\n(29/76)', ha='center', fontsize=9)
axC1.set_xticks([0,1]); axC1.set_xticklabels(['CD8-low/IDO1-low','Others'], fontsize=9.5)
axC1.set_ylabel('pCR rate (%)'); axC1.set_ylim(0,58)
axC1.set_title('C  GSE25066 (neoadjuvant chemotherapy): pCR', loc='left', fontsize=11.5, fontweight='bold')
axC1.text(-0.32, 46, "Fisher's p = 0.215 (ns)", fontsize=9, style='italic')
axC2 = axC1.inset_axes([0.46, 0.20, 0.50, 0.70])
for g,c,l in [(True,'#E41A1C','profile'), (False,'#377EB8','others')]:
    m = d66.trap.values==g
    km(d66.drfs_time.values[m].astype(float), d66.drfs_event.values[m].astype(int), axC2, c, l)
axC2.set_title('DRFS (n=117)', fontsize=8.5); axC2.set_ylim(0,1.03)
axC2.tick_params(labelsize=7.5)
axC2.text(0.97, 0.72, 'p = 0.608 (ns)\nHR 1.19 (0.61-2.30)', transform=axC2.transAxes,
          ha='right', fontsize=7.5, style='italic',
          bbox=dict(boxstyle='round', fc='white', ec='lightgrey', alpha=0.9))

axD = fig.add_subplot(gs[1,1])
x=[0,1,3,4]; vals=[20,75,17,21]
axD.bar(x, vals, color=['#E41A1C','#377EB8','#E41A1C','#377EB8'], width=0.6)
axD.text(0,23,'20%\n(1/5)', ha='center', fontsize=8.5); axD.text(1,78,'75%\n(18/24)', ha='center', fontsize=8.5)
axD.text(3,20,'17%\n(4/24)', ha='center', fontsize=8.5); axD.text(4,24,'21%\n(13/61)', ha='center', fontsize=8.5)
axD.set_xticks([0.5,3.5]); axD.set_xticklabels(['Pembrolizumab arm\n(n = 29 TNBC)','Paclitaxel control arm\n(n = 85 TNBC)'], fontsize=9)
axD.set_ylabel('pCR rate (%)'); axD.set_ylim(0,97)
axD.set_title('D  I-SPY2 (GSE194040): pCR by treatment arm', loc='left', fontsize=11.5, fontweight='bold')
axD.text(0.5, 60, 'p = 0.036\n(q = 0.083, n = 5 profile)', ha='center', fontsize=8.5)
axD.text(3.5, 30, 'p = 0.66 (ns)', ha='center', fontsize=9, style='italic')

fig.savefig(OUT, dpi=300, bbox_inches='tight')
print('saved:', OUT)
