"""Write constant/modes_HOS_SWENSE_clean.dat from constant/modes_HOS_SWENSE.dat.

For every saved record and every standard mode with dimensional k_x > KCUT
(k_x = i*pi/xlen/h, KCUT = 7.71 1/m), modeX, modeY, modeZ, modet, modeFS and
modeFSt are set to zero. For kept modes i >= 1 the potential-compatible
amplitude A = (-modeX/k + modeZ/(k tanh k))/2 (k non-dimensional) is written
back as modeX = -k A and modeZ = k tanh(k) A. Mode 0, modeY-free records,
modeAdd and modeAddt are unchanged. The output has the same size and layout.

Usage: python3 clean_modes_generator.py SRC_modes DST_modes
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hos_nwt_reader import HOSNWTModes

KCUT = 7.71
W = 18


def main(src, dst):
    f = HOSNWTModes(src)
    raw = bytearray(f.raw.encode())
    n1, n3, n2 = f.n1, f.n3_add, f.n2
    recl = W*n1
    kstar = np.arange(n1)*np.pi/f.xlen
    drop = kstar/f.h > KCUT
    keep = (~drop) & (np.arange(n1) >= 1)
    k = kstar[keep]

    def get(base):
        return np.array([float(raw[base + W*i:base + W*i + 17]) for i in range(n1)])

    def put(base, vals):
        for i in range(n1):
            raw[base + W*i:base + W*i + W] = ('%17.10E ' % vals[i]).encode()

    for blk in range(1, f.nblocks):
        base0 = blk*8*recl*n2
        rec = {nm: base0 + recl*j for j, nm in enumerate(
            ["modeX", "modeY", "modeZ", "modet", "modeFS", "modeFSt"])}
        X = get(rec["modeX"]); Z = get(rec["modeZ"])
        Xn = X.copy(); Zn = Z.copy()
        A = (-X[keep]/k + Z[keep]/(k*np.tanh(k)))/2
        Xn[keep] = -k*A
        Zn[keep] = k*np.tanh(k)*A
        Xn[drop] = 0.0; Zn[drop] = 0.0
        put(rec["modeX"], Xn)
        put(rec["modeY"], np.where(drop, 0.0, get(rec["modeY"])))
        put(rec["modeZ"], Zn)
        put(rec["modet"], np.where(drop, 0.0, get(rec["modet"])))
        put(rec["modeFS"], np.where(drop, 0.0, get(rec["modeFS"])))
        put(rec["modeFSt"], np.where(drop, 0.0, get(rec["modeFSt"])))
    with open(dst, "wb") as out:
        out.write(bytes(raw))
    print("wrote", dst, "bytes", os.path.getsize(dst), "source bytes", os.path.getsize(src))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
