# Colormodels module does not work properly on python 3.6. Use 2.7.

from colorpy import colormodels
from math import cos, sin, pi, ceil
import csv

sti_cols = []
colors = []
rgb_colors = []

for i in xrange(180):
    L = 70.0
    a = cos(i*pi/90)*40 + 5
    b = sin(i*pi/90)*40

    xyz = colormodels.xyz_from_lab([L, a, b])
    rgb = colormodels.rgb_from_xyz(xyz)
    irgb = colormodels.irgb_from_rgb(rgb)

    # convert to psychopy range (-1 - 1)
    rgb_colors.append(irgb)
    colors.append(rgb*2-1)
    if i%15 == 0:
        sti_cols.append(rgb*2-1)

# Pick X colors at max distance (including 2 slightly jittered version of that color)
no_colors = 8  # Number of colors
jitter = 12  # how much jitter

exp_colors = []
for i in range(0, 180, int(ceil(180.0/no_colors))):
    # Left probe
    if i - jitter < 0:
        exp_colors.append(rgb_colors[180 - abs(i - jitter)])
    else:
        exp_colors.append(rgb_colors[i - jitter])
    # Middle probe
    exp_colors.append(rgb_colors[i])
    # Right probe
    if i + jitter > 180:
        exp_colors.append(rgb_colors[i + jitter - 180])
    else:
        exp_colors.append(rgb_colors[i + jitter])

print(exp_colors, len(exp_colors))
with open('colors.csv', 'wb') as output:
    wr = csv.writer(output, lineterminator='\n')
    wr.writerows(exp_colors)
