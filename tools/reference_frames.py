"""Extract the two displayed reference panels from Holland et al., Supplementary Video 1."""
import argparse
import hashlib
import json
import subprocess
import tempfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = "https://pmc-oa-opendata.s3.amazonaws.com/PMC13167466.1/41592_2026_3062_MOESM3_ESM.mov"
CROPS = {"JF630b": [266, 630, 294, 328], "JF635b": [1549, 231, 294, 328]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ffmpeg", default="ffmpeg")
    parser.add_argument("--movie", type=Path, help="Previously downloaded source movie")
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as temp_dir:
        source = args.movie or Path(temp_dir) / "movie.mov"
        if not args.movie:
            urllib.request.urlretrieve(SOURCE, source)
        meta = {"source": SOURCE, "doi": "10.1038/s41592-026-03062-5",
                "credit": "Holland et al., Supplementary Video 1 (CC BY 4.0)",
                "frame": 0, "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                "transformation": "Spatial crop of the decoded first frame; publisher display intensities retained.",
                "crops": {}}
        for dye, (x, y, w, h) in CROPS.items():
            name = "h2b-" + dye.lower() + ".png"
            subprocess.run([args.ffmpeg, "-v", "error", "-y", "-i", str(source),
                            "-vf", f"format=rgb24,crop={w}:{h}:{x}:{y}", "-frames:v", "1",
                            str(ROOT / "docs" / "media" / name)], check=True)
            meta["crops"][dye] = {"image": "media/" + name, "xywh": [x, y, w, h]}
        (ROOT / "docs/data/reference-frames.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
