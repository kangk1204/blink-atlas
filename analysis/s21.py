"""Rebuild the public H2B–HaloTag, JF630b pre-illumination example.

Zenodo stores multi-gigabyte TIFFs in ZIP64. Only the first 32 MiB of the
uncompressed first movie, and the complete matching localization table and
YAML, are required here. Every whole member is checked against the ZIP CRC.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import logging
import struct
import time
import zipfile
import zlib
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

RECORD_URL = "https://zenodo.org/api/records/16941576"
DOI = "10.5281/zenodo.16941576"
PREFIX_BYTES = 32 * 1024**2
CHUNK_BYTES = 4 * 1024**2
ARCHIVE_MD5 = "500a2ce1df2c388228c150db0559edd5"
SELECTED_SHA256 = {
    "tif": "ef88d818b5c20a4d9c147e937d7c6bdae5160cd8663325abf6fe0c437aa374fd",
    "hdf5": "d9839992d9c3e0fc3203db7394ccbc5fc7afccc918ac34a3f1cd898ec1eac0ae",
    "yaml": "7de1f8c7af5dd333a560ca6dab572e89e5188826a59087cd6a332f0657ffec98",
}


def digest(blob):
    return hashlib.sha256(blob).hexdigest()


def request_bytes(url, start=None, stop=None, attempts=4):
    """Get an exact half-open byte range. Reject hosts that ignore Range."""
    headers = {"User-Agent": "H2B-Pilot-Reproducer/1.0"}
    if start is not None:
        if stop <= start:
            return b""
        headers["Range"] = f"bytes={start}-{stop - 1}"
    for attempt in range(attempts):
        try:
            with urlopen(Request(url, headers=headers), timeout=60) as response:
                if start is not None:
                    expected = f"bytes {start}-{stop - 1}/"
                    if response.status != 206 or not response.headers.get("Content-Range", "").startswith(expected):
                        raise ValueError("Zenodo did not honor the requested byte range")
                data = response.read()
                if start is not None and len(data) != stop - start:
                    raise OSError("Truncated HTTP range")
                return data
        except (HTTPError, URLError, TimeoutError, OSError) as exc:
            if attempt + 1 == attempts or isinstance(exc, HTTPError) and 400 <= exc.code < 500 and exc.code != 429:
                raise
            time.sleep(2**attempt)


class TailReader(io.RawIOBase):
    """Seekable view of the final ZIP bytes, retaining the real archive offsets."""

    def __init__(self, tail, length):
        self.tail = tail
        self.length = length
        self.pos = 0

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        self.pos = (0, self.pos, self.length)[whence] + offset
        return self.pos

    def read(self, size=-1):
        start = self.pos - (self.length - len(self.tail))
        if start < 0:
            raise OSError("ZIP directory extends beyond the fetched tail")
        data = self.tail[start:] if size < 0 else self.tail[start:start + size]
        self.pos += len(data)
        return data


def fetch(cache):
    cache.mkdir(parents=True, exist_ok=True)
    record = json.loads(request_bytes(RECORD_URL))
    archive = next(f for f in record["files"] if f["key"].startswith("Figure S21_"))
    if archive["checksum"] != "md5:" + ARCHIVE_MD5:
        raise ValueError("This S21 archive has changed since the public demo was built")
    url, length = archive["links"]["self"], archive["size"]
    tail = request_bytes(url, length - 65557, length)
    provenance = {
        "record": "https://zenodo.org/records/16941576", "doi": DOI,
        "paper": "https://doi.org/10.1002/anie.202424955",
        "title": record["metadata"]["title"],
        "creators": [creator["name"] for creator in record["metadata"]["creators"]],
        "license": record["metadata"]["license"]["id"],
        "license_url": "https://creativecommons.org/licenses/by/4.0/",
        "archive": archive["key"], "archive_bytes": length,
        "archive_md5_reported_by_zenodo": archive["checksum"].split(":")[-1],
        "retrieval": "HTTP byte ranges from the publisher's ZIP64 archive. Partial TIFF hash covers the uncompressed prefix.",
        "members": {},
    }
    with zipfile.ZipFile(TailReader(tail, length)) as z:
        members = [i for i in z.infolist() if "100ms_before_1/" in i.filename and not i.is_dir()]
    targets = {".tif": "s21-before-prefix.tif", ".hdf5": "s21-before-locs.hdf5", ".yaml": "s21-before-locs.yaml"}
    for member in members:
        suffix = Path(member.filename).suffix
        if suffix not in targets:
            continue
        header = request_bytes(url, member.header_offset, member.header_offset + 30)
        if header[:4] != b"PK\x03\x04":
            raise ValueError("Bad local ZIP header")
        name_length, extra_length = struct.unpack_from("<HH", header, 26)
        payload = member.header_offset + 30 + name_length + extra_length
        prefix = suffix == ".tif"
        wanted = min(PREFIX_BYTES, member.file_size) if prefix else member.file_size
        if member.compress_type not in (zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED):
            raise ValueError(f"Unsupported ZIP compression {member.compress_type}")
        decoder = zlib.decompressobj(-15) if member.compress_type == zipfile.ZIP_DEFLATED else None
        with io.BytesIO() as output:
            pos = 0
            while output.tell() < wanted:
                stop = min(pos + CHUNK_BYTES, member.compress_size)
                if stop <= pos:
                    raise ValueError("Truncated ZIP member")
                blob = request_bytes(url, payload + pos, payload + stop)
                chunk = decoder.decompress(blob, wanted - output.tell()) if decoder else blob[:wanted - output.tell()]
                output.write(chunk)
                pos = stop
            data = output.getvalue()
        if not prefix and zlib.crc32(data) != member.CRC:
            raise ValueError("CRC mismatch for " + member.filename)
        if digest(data) != SELECTED_SHA256[suffix[1:]]:
            raise ValueError("Source data hash mismatch for " + member.filename)
        dest = cache / targets[suffix]
        dest.write_bytes(data)
        provenance["members"][suffix[1:]] = {
            "archive_member": member.filename, "full_uncompressed_bytes": member.file_size,
            "archive_member_crc32": f"{member.CRC:08x}", "retrieved_uncompressed_bytes": len(data),
            "sha256": digest(data), "complete": len(data) == member.file_size,
        }
        print(f"{dest.name}  {len(data):,} bytes  {digest(data)}", flush=True)
    if set(provenance["members"]) != {"hdf5", "tif", "yaml"}:
        raise ValueError("The S21 archive is missing an expected member")
    (cache / "provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    return provenance


def densest_crop(table, width, height, side):
    """Choose a preview ROI solely from localizations in the 200 pilot frames."""
    import numpy as np

    good = table[(table["frame"] < 200) & (table["photons"] >= 0)]
    x = np.floor(good["x"] + 0.5).astype(int)
    y = np.floor(good["y"] + 0.5).astype(int)
    grid = np.zeros((height, width), dtype=np.int32)
    np.add.at(grid, (y, x), 1)
    integral = np.pad(grid, ((1, 0), (1, 0))).cumsum(0).cumsum(1)
    sums = integral[side:, side:] - integral[:-side, side:] - integral[side:, :-side] + integral[:-side, :-side]
    i, j = np.unravel_index(sums.argmax(), sums.shape)
    return int(j), int(i), side, side


def preview(cache, media):
    """Read raw frame strips from the prefix and overlay matching author fits."""
    import h5py
    import numpy as np
    import tifffile
    from PIL import Image, ImageDraw, ImageFont
    from analysis.pilot import picasso_calibration

    media.mkdir(parents=True, exist_ok=True)
    with h5py.File(cache / "s21-before-locs.hdf5") as hf:
        table = hf["locs"][:]
    calib = picasso_calibration(cache / "s21-before-locs.yaml")
    pixel_nm = calib["pixel_size_nm"]
    # The prefix contains complete strips and Micro-Manager frame tags. The
    # optional OME/ImageJ description IFDs live at the end of the full TIFF.
    logger = logging.getLogger("tifffile")

    def is_in_prefix(record):
        msg = record.getMessage()
        return not ("invalid value offset" in msg and any(f"TiffTag {n} @" in msg for n in (270,50838,50839)))

    logger.addFilter(is_in_prefix)
    try:
        with tifffile.TiffFile(cache / "s21-before-prefix.tif", is_ome=False, is_micromanager=False) as tif:
            frames, times = [], []
            for i in range(40):
                page = tif.pages[i]
                tags = page.tags[51123].value
                if int(tags["FrameIndex"]) != i:
                    raise ValueError("FRAME order disagrees with the Micro-Manager metadata")
                frames.append(page.asarray())
                times.append(float(tags["ElapsedTime-ms"]))
    finally:
        logger.removeFilter(is_in_prefix)
    frames = np.stack(frames)
    height, width = frames.shape[1:]
    if (height, width) != (int(calib["height_px"]), int(calib["width_px"])):
        raise ValueError("TIFF shape disagrees with Picasso's YAML")
    field_limits = np.percentile(frames[:, ::6, ::6], [3, 99.7])
    x0, y0, w, h = densest_crop(table, width, height, 224)
    crop = frames[:, y0:y0+h, x0:x0+w]
    crop_limits = np.percentile(crop, [3, 99.7])

    def grayscale(a, limits):
        values = (a.astype(float) - limits[0]) * (255. / (limits[1] - limits[0]))
        return Image.fromarray(np.clip(values, 0, 255).astype(np.uint8))

    grayscale(frames[0], field_limits).save(media / "s21-field.webp", lossless=True)
    try:
        font = ImageFont.truetype("DejaVuSans.ttf", 13)
    except OSError:
        font = ImageFont.load_default()
    views = []
    scale = 2
    for i, raw in enumerate(crop):
        im = grayscale(raw, crop_limits).resize((w*scale, h*scale), Image.Resampling.BILINEAR).convert("RGB")
        draw = ImageDraw.Draw(im)
        points = table[(table["frame"] == i) & (table["photons"] >= 0)]
        for row in points:
            x = (row["x"] + 0.5 - x0) * scale - 0.5
            y = (row["y"] + 0.5 - y0) * scale - 0.5
            if 0 <= x < im.width and 0 <= y < im.height:
                r = 5
                draw.ellipse((x-r, y-r, x+r, y+r), outline=(84, 245, 207), width=1)
        length = round(2000 / pixel_nm * scale)
        draw.rectangle((9, im.height-35, max(length+21, 68), im.height-6), fill=(14, 22, 24))
        draw.text((14, im.height-33), "2 µm", font=font, fill="white")
        draw.line((14, im.height-12, 14+length-1, im.height-12), fill="white", width=3)
        views.append(im)
    views[0].save(media / "s21-poster.webp", quality=88, method=5)
    interval = (times[-1] - times[0]) / (len(times) - 1)
    views[0].save(media / "s21-preview.webp", save_all=True, append_images=views[1:],
                  duration=round(interval), loop=0, quality=78, method=4)
    info = {
        "preview": "media/s21-preview.webp", "poster": "media/s21-poster.webp", "field": "media/s21-field.webp",
        "frame_start": 0, "frames": len(frames), "frame_interval_ms": interval,
        "interval_source": "Mean elapsed time step from Micro-Manager TIFF tag 51123 in frames 1–40",
        "crop_px": [x0, y0, w, h], "full_width_px": width, "full_height_px": height,
        "preview_width_px": views[0].width, "preview_height_px": views[0].height,
        "scale_bar_um": 2, "field_display_limits_adu": field_limits.tolist(),
        "clip_display_limits_adu": crop_limits.tolist(),
        "description": "First 40 raw frames, pre-illumination. The 224-pixel crop maximizes pilot localization count. Cyan rings mark author-fitted positions with nonnegative photons. Linear, fixed 3rd to 99.7th percentile contrast. Full-field image shows the first raw frame.",
        "attribution": "F. Walterspiel et al. (2025), Zenodo 16941576, Figure S21. CC BY 4.0. Cropped, contrast-scaled and annotated.",
    }
    return info


def render(cache, docs):
    from analysis.pilot import summarize

    source = json.loads((cache / "provenance.json").read_text(encoding="utf-8"))
    for suffix, local in [("tif", "s21-before-prefix.tif"), ("hdf5", "s21-before-locs.hdf5"), ("yaml", "s21-before-locs.yaml")]:
        actual = digest((cache / local).read_bytes())
        if actual != SELECTED_SHA256[suffix] or actual != source["members"][suffix]["sha256"]:
            raise ValueError("Cached file hash mismatch for " + local)
    source["license_url"] = "https://creativecommons.org/licenses/by/4.0/"
    movie = preview(cache, docs / "media")
    settings = {"sample_state": "live", "cell_line": "U2OS", "target": "H2B", "tag": "HaloTag",
                "reference_dye": "JF630b", "excitation_nm": 640, "frame_interval_ms": movie["frame_interval_ms"]}
    result = summarize(cache / "s21-before-locs.hdf5", settings, min_photons=0)
    result.update(source=source, media=movie, intended_use="Live-cell measurement and forecast QC, outside the fixed COS-7 dye-ranking protocol.")
    (docs / "data").mkdir(parents=True, exist_ok=True)
    (docs / "data" / "pilot.json").write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({k: result["evaluation"][k] for k in ["model_mae", "baseline_mae", "model_improved"]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--docs", type=Path, default=Path("docs"))
    parser.add_argument("--reuse", action="store_true", help="Check locally cached source hashes and skip fetching")
    parser.add_argument("--fetch-only", action="store_true", help="Fetch upstream inputs without rendering")
    args = parser.parse_args()
    if not args.reuse:
        fetch(args.cache)
    if not args.fetch_only:
        render(args.cache, args.docs)
