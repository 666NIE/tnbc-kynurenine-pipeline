# tools/make_figS6_concept.py -- 图S6 双模型概念图（matplotlib 手绘示意图，非数据驱动）
# 用法: python tools/make_figS6_concept.py   （仅需 matplotlib；产物输出到 figures/）
# 注意: 本图为示意性机制图，数值均取自 manuscript 终值；非 run_all 管线组成部分
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "figures")

fig, ax = plt.subplots(figsize=(16, 10), dpi=200)
ax.set_xlim(0, 160); ax.set_ylim(0, 100); ax.axis("off")

C_BLUE_D, C_BLUE_L = "#2E5E8C", "#E8F0F8"
C_OR_D,  C_OR_L  = "#B05A1E", "#FBEEE2"
C_HUB_D, C_HUB_L = "#7B1F1F", "#F7E3E3"
C_GREY_D = "#4A4A4A"
C_GRN_D = "#2E7D4F"

def rbox(x, y, w, h, fc, ec, lw=1.6, r=1.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0.02,rounding_size={r}",
                                fc=fc, ec=ec, lw=lw, zorder=2))
def header(x, y, w, h, fc, text, fs=12.5):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=1.0",
                                fc=fc, ec="none", zorder=3))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs,
            fontweight="bold", color="white", zorder=4)
def bullets(x, y, items, color, fs=9.8, dy=5.4, msize=7):
    for i, it in enumerate(items):
        yy = y - i * dy
        ax.plot([x + 1.2], [yy], marker="o", ms=msize, color=color, zorder=4)
        ax.text(x + 3.4, yy, it, ha="left", va="center", fontsize=fs, color="#222222", zorder=4)
def strip(x, y, w, h, fc, text, fs=9.5, tc=None):
    ax.add_patch(Rectangle((x, y), w, h, fc=fc, ec="none", zorder=3))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs,
            fontweight="bold", color=tc or "#222222", zorder=4)

ax.text(80, 96, "Kynurenine-Centered Metabolic\u2013Immune Trap: the Dual-Arm Model",
        ha="center", fontsize=19, fontweight="bold", color="#1a1a1a")
ax.text(80, 91.8, "TCGA discovery (n = 123)  \u2192  METABRIC validation (n = 320)  \u2192  "
        "single-cell & spatial resolution  \u2192  concordant prediction frameworks",
        ha="center", fontsize=11, color="#666666")

rbox(6, 38, 48, 46, "white", C_BLUE_D, lw=2)
header(6, 76.5, 48, 7.5, C_BLUE_D, "UPSTREAM INFLAMMATORY ARM  (IFN-\u03b3\u2013induced)")
bullets(9, 72.5, [
    "IFN-\u03b3 drives IDO1 / TDO2 transcriptional activation",
    "CD8\u207a T-cell-rich, inflamed microenvironment",
    "IDO1 high  \u2192  favorable OS (HR 0.78, p = 0.0018)",
    "Checkpoint-blockade\u2013responsive tendency\n(TIDE / IMPRES)"], C_BLUE_D)
strip(6, 38.5, 48, 6, C_BLUE_L, "TCGA C3 (n = 50)   \u00b7   METABRIC M2+M3 (n = 223)", tc=C_BLUE_D)

rbox(106, 38, 48, 46, "white", C_OR_D, lw=2)
header(106, 76.5, 48, 7.5, C_OR_D, "DOWNSTREAM METABOLIC ARM  (constitutive)")
bullets(109, 72.5, [
    "KP program weight \u2191 (relative ssGSEA ranking)",
    "Rate-limiting enzymes not elevated; downstream\nnodes AFMID / AADAT / GPT2 show up-trends",
    "Myeloid KP flux: IL4I1 / KYNU / KMO in\nmacrophages & DCs (T cells < 1.5%)",
    "TAM M2\u2013linked exclusion signature (TIDE)"], C_OR_D)
strip(106, 38.5, 48, 6, C_OR_L, "TCGA C1 (n = 28)   \u00b7   METABRIC M1 (n = 97)   \u00b7   OS HR 1.49, p = 0.033", tc=C_OR_D)

rbox(60, 46, 40, 32, C_HUB_L, C_HUB_D, lw=2.2, r=1.6)
ax.text(80, 71.5, "KYNURERENINE-CENTERED\nMETABOLIC\u2013IMMUNE TRAP", ha="center", va="center",
        fontsize=12.5, fontweight="bold", color=C_HUB_D, zorder=4)
ax.plot([63, 97], [64.5, 64.5], color=C_HUB_D, lw=0.8, zorder=4)
ax.text(80, 57.5, "Score\u2013expression discordance\nKyn score \u2191   \u00b7   IDO1 expression \u2193\nCD8 / cytotoxicity \u2193",
        ha="center", va="center", fontsize=10, color="#3a3a3a", zorder=4)

ax.add_patch(FancyArrowPatch((54.2, 61), (59.8, 61), arrowstyle="-|>", mutation_scale=26, lw=2.6, color=C_BLUE_D, zorder=1))
ax.add_patch(FancyArrowPatch((105.8, 61), (100.2, 61), arrowstyle="-|>", mutation_scale=26, lw=2.6, color=C_OR_D, zorder=1))
ax.text(57, 63.5, "feeds", fontsize=8.5, color=C_BLUE_D, ha="center")
ax.text(103, 63.5, "feeds", fontsize=8.5, color=C_OR_D, ha="center")

rbox(6, 24, 148, 9.5, "white", C_GREY_D, lw=1.4)
ax.text(80, 30.4, "VALIDATION CHAIN", ha="center", fontsize=10, fontweight="bold", color=C_GREY_D, zorder=4)
ax.text(80, 26.6, "consensus clustering (k = 3)  \u2192  dual-track sensitivity  \u2192  swap-in clustering (original C1: 0% in abs-high cluster)"
        "  \u2192  uniform-metric ssGSEA (100% recovery)  \u2192  scRNA atlas (n = 100,064 cells)  \u2192  spatial Visium (n = 11,877 spots)  \u2192  TIDE + IMPRES concordance",
        ha="center", fontsize=8.8, color="#3a3a3a", zorder=4)

rbox(6, 8, 148, 11.5, "white", C_GRN_D, lw=1.4)
ax.text(80, 16.4, "CANDIDATE TRIPLEX IHC PANEL  (research-stage stratification)", ha="center",
        fontsize=10, fontweight="bold", color=C_GRN_D, zorder=4)
ax.text(80, 11.6, "KP downstream nodes (AFMID / AADAT / GPT2) + CD8 + IDO1      \u2192      IDO1-high: inflammation-induced type (ICB benefit)"
        "      |      IDO1-low + downstream-high + CD8-low: constitutive trap (anti-exclusion strategy)      \u2014      prospective validation pending",
        ha="center", fontsize=8.8, color="#3a3a3a", zorder=4)

plt.tight_layout()
fig.savefig(os.path.join(OUT, "FigureS6_dual_model_v2.png"), dpi=200, bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(OUT, "FigureS6_dual_model_v2.pdf"), bbox_inches="tight", facecolor="white")
print("saved to", OUT)
