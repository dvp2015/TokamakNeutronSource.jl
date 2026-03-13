# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
from typing import Generator


# %%
from pathlib import Path

import numpy as np

from tqdm import tqdm

import mckit as mc


# %%
# !pwd


# %%
def read_csv(path):
    return np.genfromtxt(path, dtype=float, delimiter=",")


# %%
rbins = read_csv("../wrk/rbins.csv") * 100.0  # m -> cm
zbins = read_csv("../wrk/zbins.csv") * 100.0
intensity = read_csv("../wrk/src.csv")


# %%
rbins.shape, zbins.shape, intensity.shape


# %%
from mckit import source
from mckit.source import create_bin_distributions


start_name=1
aux_name = start_name + 10
aux_name, x_distr = create_bin_distributions(rbins, aux_name)
aux_name, y_distr = create_bin_distributions(zbins, aux_name)
probs = []
r_indices = []
z_indices = []
Nx = len(x_distr) - 1
Ny = len(y_distr) - 1 

for i in range(Nx):
    for j in range(Ny):
        if intensity[i, j] > 0:
            r_indices.append(x_distr[i])
            z_indices.append(y_distr[j])
            probs.append(intensity[i, j])

r_dist = source.Distribution(start_name + 1, r_indices, probs, 'RAD')
z_dist = source.Distribution(start_name + 2, z_indices, r_dist, 'EXT')
src = source.Source(RAD=r_dist, EXT=z_dist)


# %%
with open("../wrk/sdef-dd-2.txt", "w") as io:    
    first_line, text = src.mcnp_repr().split("\n", maxsplit=1)
    comment_line = "c DD source from TRT_215_8T_NBI plasma model"
    first_line = first_line.lower().replace("rad=d2", "cell=47 erg=d1 axs=0. 0. 1. rad=d2")
    second_line = "sp1  -4  -0.01  -2  $ DD fusion 10 keV"
    io.write("\n".join([comment_line, first_line, second_line, text, ""]))


# %%
def ptrac_to_csv(dst: Path, src: Path) -> None:
    PREV_LINE = "       9000         1        40        47         0"
    wait_prev = True
    with src.open() as inp, dst.open("w") as out:
        with tqdm(total=1000000) as pbar:
            for line in inp:
                if wait_prev:
                    if line.startswith(PREV_LINE):
                        wait_prev = False
                else:
                    x, y, z = (float(x) for x in line.split())
                    r = np.linalg.norm([x, y])
                    print(r, z, file=out)
                    pbar.update()
                    wait_prev = True


# %%
ptrac_to_csv(Path("../wrk/ptrac.csv"), Path("../wrk/ptrac"))


# %%


datap = np.genfromtxt("../wrk/ptrac.csv")


# %%


datap

