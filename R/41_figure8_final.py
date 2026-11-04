# 41_figure8_final.py — Figure 8（跨队列验证）终版生成脚本
# 与投稿版 Figure 8 一一对应。数据来源：
#   results/36_gse21653_perpatient.csv   （GSE21653 逐患者，由 36_gse21653_decomposition.R 导出）
#   results/40_gse25066_trap_assignment.csv（GSE25066 逐患者，由 40_gse25066_validation.R 输出）
# A/D 面板数据为已发表汇总值（见正文 Figure 8 图注）。用法：python3 41_figure8_final.py
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2

# ---------- KM 工具 ----------
def km_curve(t, e):
    order = np.argsort(t); t, e = np.asarray(t)[order], np.asarray(e)[order]
    ts = np.unique(t[e == 1]); S = 1.0
    pt, ps, cx, cy = [0], [1.0], [], []
    for tt in ts:
        nr = (t >= tt).sum(); de = ((t == tt) & (e == 1)).sum()
        S *= (1 - de / nr)
        pt += [tt, tt]; ps += [ps[-1], S]
        for tc in t[(t == tt) & (e == 0)]: cx.append(tc); cy.append(S)
    pt.append(t.max() + 0.01); ps.append(S)
    return pt, ps, cx, cy

def logrank(t, e, grp):
    t = np.asarray(t); e = np.asarray(e); grp = np.asarray(grp)
    ts = np.unique(t[e == 1]); O = E = V = 0.0
    for tt in ts:
        n = (t >= tt).sum(); n1 = ((t >= tt) & grp).sum()
        d = ((t == tt) & (e == 1)).sum(); d1 = ((t == tt) & (e == 1) & grp).sum()
        if n <= 1: continue
        E += d * n1 / n; V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1); O += d1
    return 1 - chi2.cdf((O - E) ** 2 / V, 1)

# ---------- 数据 ----------
d53 = pd.read_csv("results/36_gse21653_perpatient.csv")
tr_trap = km_curve(d53.t[d53.trap], d53.e[d53.trap])
tr_oth  = km_curve(d53.t[~d53.trap], d53.e[~d53.trap])

g66 = pd.read_csv("results/40_gse25066_trap_assignment.csv").dropna(subset=["drfs_time", "drfs_event"]).query("drfs_time > 0")
c_tr = km_curve(g66.drfs_time[g66.trap], g66.drfs_event[g66.trap])
c_ot = km_curve(g66.drfs_time[~g66.trap], g66.drfs_event[~g66.trap])

# ---------- 画图 ----------
fig = plt.figure(figsize=(14, 10.5))
gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.26)
T = dict(loc='left', fontsize=11.5, fontweight='bold')

axA = fig.add_subplot(gs[0, 0])
coh = ['TCGA C1\n(n=123)', 'METABRIC M1\n(n=320)', 'GSE21653\n(n=85)', 'I-SPY2\n(n=114 TNBC)', 'GSE25066\n(n=182)']
pct = [23, 30, 28.2, 25.4, 34.6]
axA.bar(range(5), pct, color=['#4C72B0', '#4C72B0', '#55A868', '#DD8452', '#C44E52'])
for i, v in enumerate(pct): axA.text(i, v + 0.8, f'{v}%', ha='center', fontsize=9.5, fontweight='bold')
axA.axhspan(23, 35, color='grey', alpha=0.12)
axA.set_xticks(range(5)); axA.set_xticklabels(coh, fontsize=8.5)
axA.set_ylabel('CD8-low/IDO1-low profile (%)'); axA.set_ylim(0, 46)
axA.set_title('A  CD8-low/IDO1-low profile prevalence across five cohorts', **T)

axB = fig.add_subplot(gs[0, 1])
axB.plot(*tr_trap[:2], drawstyle='steps-post', color='#E41A1C', lw=1.8, label='CD8-low/IDO1-low (n=24)')
axB.plot(*tr_oth[:2], drawstyle='steps-post', color='#377EB8', lw=1.8, label='Others (n=61)')
axB.plot(tr_trap[2], tr_trap[3], '|', color='#E41A1C', ms=8, mew=1.2)
axB.plot(tr_oth[2], tr_oth[3], '|', color='#377EB8', ms=8, mew=1.2)
axB.set_xlim(-0.32, 5.4); axB.set_ylim(-0.035, 1.06)
axB.set_title('B  GSE21653: disease-free survival', **T)
axB.set_xlabel('Time (years)'); axB.set_ylabel('DFS probability')
axB.grid(True, alpha=0.3); axB.legend(fontsize=9, loc='lower left')
axB.text(0.96, 0.60, 'log-rank p < 0.001\nHR 4.84 (95% CI 2.23-10.50)\nn = 85, 27 events',
         transform=axB.transAxes, ha='right', fontsize=9,
         bbox=dict(boxstyle='round', fc='white', ec='grey', alpha=0.9))

axC1 = fig.add_subplot(gs[1, 0])
axC1.bar([0, 1], [25.6, 38.2], color=['#E41A1C', '#377EB8'], width=0.55)
axC1.text(0, 27.5, '25.6%\n(10/39)', ha='center', fontsize=9)
axC1.text(1, 40.5, '38.2%\n(29/76)', ha='center', fontsize=9)
axC1.set_xticks([0, 1]); axC1.set_xticklabels(['CD8-low/IDO1-low', 'Others'], fontsize=9.5)
axC1.set_ylabel('pCR rate (%)'); axC1.set_ylim(0, 58)
axC1.set_title('C  GSE25066 (neoadjuvant chemotherapy): pCR', **T)
axC1.text(-0.32, 47, "Fisher's p = 0.215 (ns)", fontsize=9, style='italic')

axC2 = axC1.inset_axes([0.44, 0.16, 0.52, 0.76])
axC2.plot(*c_tr[:2], drawstyle='steps-post', color='#E41A1C', lw=1.4, label='CD8-low/IDO1-low')
axC2.plot(*c_ot[:2], drawstyle='steps-post', color='#377EB8', lw=1.4, label='Others')
axC2.plot(c_tr[2], c_tr[3], '|', color='#E41A1C', ms=5, mew=1)
axC2.plot(c_ot[2], c_ot[3], '|', color='#377EB8', ms=5, mew=1)
axC2.set_title('DRFS (n = 117)', fontsize=9, fontweight='bold')
axC2.set_ylim(-0.05, 1.08); axC2.set_xlim(-0.3, 8.2)
axC2.tick_params(labelsize=7)
axC2.set_xlabel('Years', fontsize=8); axC2.set_ylabel('DRFS', fontsize=8)
axC2.grid(True, alpha=0.25); axC2.legend(fontsize=6.5, loc='lower left')
axC2.text(0.97, 0.62, 'p = 0.608 (ns)\nHR 1.19 (0.61-2.30)',
          transform=axC2.transAxes, ha='right', fontsize=7.5,
          bbox=dict(boxstyle='round', fc='white', ec='grey', alpha=0.9))

axD = fig.add_subplot(gs[1, 1])
x = [0, 1, 3, 4]; vals = [20, 75, 17, 21]
axD.bar(x, vals, color=['#E41A1C', '#377EB8', '#E41A1C', '#377EB8'], width=0.6)
axD.text(0, 23, '20%\n(1/5)', ha='center', fontsize=8.5); axD.text(1, 78, '75%\n(18/24)', ha='center', fontsize=8.5)
axD.text(3, 20, '17%\n(4/24)', ha='center', fontsize=8.5); axD.text(4, 24, '21%\n(13/61)', ha='center', fontsize=8.5)
axD.set_xticks([0.5, 3.5]); axD.set_xticklabels(['Pembrolizumab arm\n(n = 29 TNBC)', 'Paclitaxel control arm\n(n = 85 TNBC)'], fontsize=9)
axD.set_ylabel('pCR rate (%)'); axD.set_ylim(0, 97)
axD.set_title('D  I-SPY2 (GSE194040): pCR by treatment arm', **T)
axD.text(0.5, 60, 'p = 0.036\n(q = 0.083, n = 5 profile)', ha='center', fontsize=8.5)
axD.text(2.28, 52, 'p = 0.66 (ns)', ha='left', fontsize=9, style='italic')

fig.savefig("figures/Figure8.pdf", bbox_inches='tight')
fig.savefig("figures/Figure8.png", dpi=300, bbox_inches='tight')
print("done -> figures/Figure8.pdf")
