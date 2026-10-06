"""Design scan: achievable parameter accuracy (Cramer-Rao bound) vs instrument
noise level, for each calibration problem. Produces the instrument-selection
table of the experiment design report.
"""

import numpy as np
from calib_problems import PROBLEMS, theta_true
from calib_invert import fim_analysis

NOISE_LEVELS = [0.001, 0.005, 0.01, 0.02, 0.05]   # relative std

print(f'{"problem":22s} {"param":9s} ' + ' '.join(f'{n*100:8.1f}%' for n in NOISE_LEVELS))
for name, prob in PROBLEMS.items():
    base_noise = list(prob['obs_noise'])
    for j, (pn, *_rest) in enumerate(prob['params']):
        row = []
        for nl in NOISE_LEVELS:
            prob['obs_noise'] = [(g, nl) for g, _s in base_noise]
            fa = fim_analysis(prob, theta_true(prob))
            row.append(fa['rel_std'][j])
        prob['obs_noise'] = list(base_noise)
        tag = name if j == 0 else ''
        print(f'{tag:22s} {pn:9s} ' + ' '.join(f'{r*100:8.2f}%' for r in row))
