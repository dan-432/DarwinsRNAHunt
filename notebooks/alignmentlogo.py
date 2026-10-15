import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import logomaker


# # palette.py
KEA = {
    "green":       "#6C803A",
    "dark_brown":  "#7B5C34",
    "light_brown": "#AB7C47",
    "yellow":      "#CCAE42",
    "orange":      "#D73202",
    "charcoal":    "#272318",
    "light_grey":  "#D3CDBF",
}

_GROUPS = {
    "AVLIM": KEA["green"],        # hydrophobic
    "FWY":   KEA["dark_brown"],   # aromatic
    "STNQ":  KEA["light_brown"],  # polar
    "KRH":   KEA["charcoal"],     # positive
    "DE":    KEA["orange"],       # negative
    "GPC":   KEA["yellow"],       # special
}

KEA_LOGO = {aa: colour for group, colour in _GROUPS.items() for aa in group}

def read_fasta(path):
    seqs, cur = [], []
    with open(path) as fasta:
        for line in fasta:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                cur = []
            else:
                cur.append(line.strip().upper())
    if cur:
        seqs.append("".join(cur))
    return seqs

seqs = read_fasta("results/00_controls/rfam/RF00050/05_alignments/cluster_204.msa")

# seqs = list of aligned sequences, as in the earlier example
start, end = 411, 431  # 0-based, end-exclusive; pick a conserved stretch

counts = logomaker.alignment_to_matrix(seqs, to_type="counts")
info = logomaker.transform_matrix(counts, from_type="counts", to_type="information")
window = info.iloc[start:end]

fig, ax = plt.subplots(figsize=(4, 1.8))
logo = logomaker.Logo(window, ax=ax, color_scheme=KEA_LOGO)
ax.set_ylabel("bits")
ax.set_xlabel("alignment column")
ax.set_ylim(0, 4.32)  # log2(20), the maximum for a protein alphabet
fig.tight_layout()
fig.savefig("logo_window.pdf")  # vector output is better for a thesis




START, END = 411, 431  # 0-based, end-exclusive; pick a conserved stretch



def text_colour(hex_colour):
    """White text on dark cells, charcoal on light ones."""
    r, g, b = (int(hex_colour[i:i + 2], 16) / 255 for i in (1, 3, 5))
    return "white" if 0.299 * r + 0.587 * g + 0.114 * b < 0.55 else KEA["charcoal"]


def pick_sequences(seqs, start, end, n_show):
    """Gap-free rows in the window first, then evenly spaced down the alignment."""
    full = [s for s in seqs if "-" not in s[start:end]]
    pool = full if len(full) >= n_show else seqs
    step = max(len(pool) // n_show, 1)
    return pool[::step][:n_show]


def draw_alignment(ax, window_seqs):
    for row, seq in enumerate(window_seqs):
        for col, aa in enumerate(seq):
            if aa == "-":
                ax.text(col, row, "-", ha="center", va="center",
                        color=KEA["light_grey"], fontsize=7)
                continue
            colour = KEA_LOGO.get(aa, KEA["light_grey"])
            ax.add_patch(Rectangle((col - 0.5, row - 0.5), 1, 1,
                                   facecolor=colour, edgecolor="white", linewidth=0.8))
            ax.text(col, row, aa, ha="center", va="center", fontsize=7,
                    family="monospace", color=text_colour(colour))
    ax.set_ylim(len(window_seqs) - 0.5, -0.5)
    ax.set_yticks([])
    ax.tick_params(bottom=False, labelbottom=False)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)


def alignment_and_logo(seqs, start, end, n_show=8, out="alignment_logo.pdf"):
    shown = pick_sequences(seqs, start, end, n_show)
    window = [s[start:end] for s in shown]
    width = end - start

    counts = logomaker.alignment_to_matrix(window, to_type="counts")
    info = logomaker.transform_matrix(counts, from_type="counts", to_type="information")

    fig, (ax_aln, ax_logo) = plt.subplots(
        2, 1, figsize=(width * 0.28, n_show * 0.2 + 1.8),
        gridspec_kw={"height_ratios": [n_show, 5], "hspace": 0.12}, sharex=True,
    )
    draw_alignment(ax_aln, window)

    counts = logomaker.alignment_to_matrix(seqs, to_type="counts")
    info = logomaker.transform_matrix(counts, from_type="counts", to_type="information")
    window = info.iloc[start:end]

    logomaker.Logo(info, ax=ax_logo, color_scheme=KEA_LOGO)
    ax_logo.set_ylim(0, 2)
    ax_logo.set_ylabel("bits")
    ticks = list(range(0, width, 5))
    ax_logo.set_xticks(ticks)
    ax_logo.set_xticklabels([str(start + t + 1) for t in ticks])  # 1-based alignment columns
    ax_logo.set_xlabel("alignment column")
    ax_logo.set_xlim(-0.5, width - 0.5)

    fig.savefig(out, bbox_inches="tight")


seqs = read_fasta("results/00_controls/rfam/RF00050/05_alignments/cluster_204.msa")  # from the first example
alignment_and_logo(seqs, start=START, end=END)  # use the window from the earlier script