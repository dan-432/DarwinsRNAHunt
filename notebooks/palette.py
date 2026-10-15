# palette.py
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