"""Checks of an installation that its tests on the source cannot see: the built app leaves out parts of its libraries
(packaging/stannotator.spec), and a new version of a library can move what the app needs into a part it leaves out
(imagecodecs 2026.11.11 renamed its JPEG module, and the app could neither read nor write JPEG tiles).

run() writes and reads back a small tiled TIFF in every compression of the H&E images the app opens (whole-slide
OME-TIFF, SVS, NDPI) and of the image pyramids it makes, through tifffile and imagecodecs as the app does. The
release builds run it on the built app (stagent selfcheck) and fail when a compression does not work."""
import io
import sys

import numpy as np

# tifffile's names of the compressions: JPEG (the app's pyramids, most whole-slide images), JPEG 2000 (SVS), and the
# lossless ones of OME-TIFF exports
COMPRESSIONS = ["jpeg", "jpeg2000", "zlib", "lzw", "zstd", "webp", "png"]


def run(out=print):
    import imagecodecs
    import tifffile
    y, x = np.mgrid[:256, :256]
    img = np.stack([180 + x // 8, 150 + y // 8, 200 - (x + y) // 16], -1).astype(np.uint8)  # smooth (JPEG is lossy)
    out(f"tifffile {tifffile.__version__}, imagecodecs {imagecodecs.__version__}")
    bad = []
    for c in COMPRESSIONS:
        try:
            b = io.BytesIO()
            tifffile.imwrite(b, img, tile=(128, 128), photometric="rgb", compression=c)
            b.seek(0)
            back = tifffile.imread(b)
            err = float(np.abs(back.astype(float) - img).mean()) if back.shape == img.shape else np.inf
            if err > 3:
                raise ValueError(f"read back differs (shape {back.shape}, mean difference {err:.1f})")
            out(f"  {c}: ok")
        except Exception as e:
            out(f"  {c}: FAILED: {type(e).__name__}: {e}")
            bad.append(c)
    if bad:
        out("image compressions that do not work: " + ", ".join(bad))
        return 1
    out("selfcheck ok")
    return 0


if __name__ == "__main__":
    sys.exit(run())
