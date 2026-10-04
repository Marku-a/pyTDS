# Reel audio credits

## Voice-over
- Engine: Piper TTS (piper-tts 1.8.0 from PyPI, OHF-Voice/rhasspy), run offline.
- Voice: `en-us-libritts-high` (rhasspy/piper release v0.0.2), speaker id 6555, length_scale 1.3, noise_scale 0.5, noise_w 0.6, 0.45 s of silence between sentences.
- Dataset: LibriTTS (http://www.openslr.org/60/), licence CC BY 4.0 (public and commercial use allowed with attribution). Per the voice's MODEL_CARD.
- Attribution: "Voice: Piper en-us-libritts-high (speaker 6555), trained on LibriTTS (Zen et al., 2019), CC BY 4.0, https://www.openslr.org/60/"
- Selection: 61 speakers (random sample, seed 0, plus the former speaker 4535) ranked by word error rate of a CMU pocketsphinx (en-us) recogniser on 3 script sentences; speaker 6555 was 2nd by WER (0.13; former speaker 4535: 0.41), had near-lowest spectral flatness and a speaking rate of about 2.2-2.4 words/s.
- Processing (ffmpeg): highpass 80 Hz, +2.5 dB presence at 3 kHz, light compression (ratio 3, -20 dB), gain to -16 LUFS, peak limiter.

Candidates rejected on licence: Ryan (CC BY-NC-SA 4.0), Lessac (Blizzard 2013, research licence), Amy/Alan (MODEL_CARD says only "See URL", not verified).
Kathleen (CC0, low quality, 16 kHz) was evaluated as a candidate (WER 0.17 on the same test) but not chosen.

## Music
Original, synthesised with numpy by `make_music.py`; no third-party material.
