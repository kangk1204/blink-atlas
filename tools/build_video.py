"""Mux a 120-second browser recording with English narration and timed captions."""
import argparse
import json
import subprocess
import tempfile
import wave
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RATE = 24000


def stamp(seconds):
    ms = round(seconds * 1000)
    return f"{ms//3600000:02d}:{ms//60000%60:02d}:{ms//1000%60:02d}.{ms%1000:03d}"


def duration(path):
    with wave.open(str(path), "rb") as inp:
        return inp.getnframes() / inp.getframerate()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("capture", type=Path)
    parser.add_argument("wav_dir", type=Path)
    parser.add_argument("--ffmpeg", default="ffmpeg")
    args = parser.parse_args()
    media = ROOT / "docs/media"
    scenes = json.loads((ROOT / "tools/tour.json").read_text(encoding="utf-8"))
    captions = ["WEBVTT", ""]
    length = scenes[-1]["end"]
    pcm = bytearray(round(RATE * length) * 2)
    for i, scene in enumerate(scenes):
        files = [args.wav_dir / f"{i:02d}-{j:02d}.wav" for j in range(len(scene["cues"]))]
        available = scene["end"] - scene["start"] - 0.8 - 0.22 * (len(files)-1)
        speed = max(0.92, sum(duration(f) for f in files) / available)
        offset = round((scene["start"] + 0.3) * RATE)
        for path, cue in zip(files, scene["cues"]):
            audio = subprocess.check_output([args.ffmpeg, "-v", "error", "-i", str(path),
                      "-af", f"atempo={speed:.6f}", "-ac","1","-ar",str(RATE),"-f","s16le","pipe:1"])
            end = offset + len(audio)//2
            if end > scene["end"] * RATE:
                raise ValueError("Voiceover exceeds the scene boundary")
            pcm[offset*2:end*2] = audio
            captions.extend([f"{stamp(offset/RATE)} --> {stamp(end/RATE)}", cue, ""])
            offset = end + round(0.22 * RATE)
    vtt = media / "walkthrough.vtt"
    vtt.write_text("\n".join(captions).rstrip() + "\n", encoding="utf-8")
    with tempfile.TemporaryDirectory() as folder:
        voice = Path(folder) / "voice.wav"
        with wave.open(str(voice), "wb") as out:
            out.setparams((1, 2, RATE, 0, "NONE", "not compressed"))
            out.writeframes(pcm)
        subprocess.run([args.ffmpeg, "-v","error","-y","-ss","0.35","-i",str(args.capture),
                      "-i",str(voice),"-i",str(vtt),"-map","0:v:0","-map","1:a:0","-map","2:s:0",
                      "-t",str(length),"-vf","fps=24,format=yuv420p","-c:v","libx264","-preset","fast",
                      "-crf","20","-c:a","aac","-b:a","96k","-c:s","mov_text",
                      "-metadata:s:s:0","language=eng","-metadata:s:a:0","language=eng",
                      "-metadata","title=Blink Atlas walkthrough","-movflags","+faststart",
                      str(media / "walkthrough.mp4")],check=True)


if __name__ == "__main__":
    main()
