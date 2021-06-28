from ._tamis.ebcd import ebcd_core


def ebcd(signal, weights, lbd):
    return ebcd_core(signal, weights, lbd)
