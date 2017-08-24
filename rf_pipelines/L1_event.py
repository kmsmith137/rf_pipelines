import sys
import numpy as np


def get_dtypes(src):
    """
    """
    config = {}
    if isinstance(src, str):
        execfile(src, config)
        config.setdefault('ntrees', 1)
        config.setdefault('nbeta', [1])
        config.setdefault('nsm', [1])
        for name in 'nbeta', 'nsm':
            if np.isscalar(config[name]):
                config[name] = [config[name]]
    else:
        config = src.config

    compact = [('beam', np.uint16),
               ('itree', np.uint8),
               ('snr', np.float32),
               ('time', 'datetime64[us]'),
               ('dm', np.float32),
               ('snr_vs_dm', np.float32, 17),
               ('grade', np.uint8)]

    if config['ntrees'] > 1:
        compact.append(('snr_vs_itree', np.float32, config['ntrees']))

    if max(config['nbeta']) > 1:
        compact.append(('beta', np.uint8))
        compact.append(('snr_vs_beta', np.float32, max(config['nbeta'])))

    if max(config['nsm']) > 1:
        compact.append(('sm', np.uint8))
        compact.append(('snr_vs_sm', np.float32, max(config['nsm'])))

    extended = compact + [('rajd', np.float32),
                          ('decjd', np.float32)]

    return compact, extended

try:
    from chime_frb_L2_L3_config import L1_bonsai_config
    L1_HDR_DTYPE, L1_EXTENDED = get_dtypes(L1_bonsai_config)
except:
    pass
