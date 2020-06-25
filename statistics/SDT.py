import math
from scipy.stats import norm


def SDT(hits, misses, fas, crs):
    """ To calculate the most important Signal detection theory values d', beta, c and Ad',
    using percentile point function (ppf from scipy).

    Parameters
    ----------
    hits : [type]
        hit
    misses : [type]
        misses
    fas : [type]
        false alarms
    crs : [type]
        correct rejections

    Returns
    -------
    out : dict
        Contains d (d-prime: sensitivity), beta (observer criterion), C (observer criterion), 
        and Ad (A d-prime: sensitivity for present/absent task)
    """

    # Floors an ceilings are replaced by half hits and half FA's
    half_hit = 0.5 / (hits + misses)
    half_fa = 0.5 / (fas + crs)

    # Calculate hit_rate and avoid d' infinity
    hit_rate = hits / (hits + misses)
    if hit_rate == 1:
        hit_rate = 1 - half_hit
    if hit_rate == 0:
        hit_rate = half_hit

    # Calculate false alarm rate and avoid d' infinity
    fa_rate = fas / (fas + crs)
    if fa_rate == 1:
        fa_rate = 1 - half_fa
    if fa_rate == 0:
        fa_rate = half_fa

    # Return d', beta, c and Ad'
    out = {}
    Z = norm.ppf
    out['d'] = Z(hit_rate) - Z(fa_rate)
    out['beta'] = math.exp((Z(fa_rate)**2 - Z(hit_rate)**2) / 2)
    out['c'] = -(Z(hit_rate) + Z(fa_rate)) / 2
    out['Ad'] = norm.cdf(out['d'] / math.sqrt(2))

    return out
