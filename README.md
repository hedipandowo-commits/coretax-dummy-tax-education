# Data dan Skrip Analisis — Sistem Dummy Coretax sebagai Media Pembelajaran Administrasi Perpajakan

Repositori pendukung artikel *"Perancangan dan Penerapan Sistem Dummy Coretax sebagai Media Pembelajaran Administrasi Perpajakan: Uji Efektivitas, Kebergunaan, dan Penerimaan Teknologi"* (Pandowo, Sugiharto, Amrullah; Politeknik Negeri Madiun).

Isi repositori memungkinkan seluruh angka inferensial pada artikel (H1–H14) direproduksi.

## Struktur
```
data/       data_respons_individual.csv (N = 96, anonim) + KAMUS_DATA.md
code/       analisis_statistik.py   -> uji-t, Wilcoxon, SUS, PLS-SEM, bootstrap, mediasi, moderasi
            buat_gambar.py          -> Gambar 2 dan 3 dari hasil analisis
notebook/   Analisis_Coretax_Colab.ipynb -> menjalankan semuanya di Google Colab
results/    hasil_analisis.json     -> keluaran lengkap (angka pada Tabel 3–9)
figures/    gambar2_model_tax_tam.png, gambar3_simple_slopes.png
```

## Cara menjalankan
**Google Colab (tanpa instalasi):** buka `notebook/Analisis_Coretax_Colab.ipynb` di Colab, lalu *Runtime → Run all*. Waktu proses sekitar 1–3 menit.

**Lokal:**
```bash
pip install -r requirements.txt
python code/analisis_statistik.py        # menulis results/hasil_analisis.json
python code/buat_gambar.py               # menulis figures/*.png
```
Opsi: `--boot 1000` untuk uji cepat (artikel memakai 5.000 sampel ulang), `--data` dan `--out` untuk jalur kustom.

## Metode (ringkas)
- H1: uji-t berpasangan dan Wilcoxon (satu sisi), dz, N-Gain.
- H2, H3a, H3b: uji-t satu sampel dan Wilcoxon terhadap 68 (SUS) dan 3,0 (PU, PEOU).
- H4–H14: PLS-SEM mode A, skema pembobotan jalur, data terstandardisasi, konvergensi 1e-7, bootstrap 5.000 (persentil 95%), mediasi (Zhao et al., 2010), moderasi dua tahap.
- Seed acak tetap (`20260110`), sehingga hasil bootstrap identik pada setiap eksekusi dengan versi pustaka yang sama.
- Implementasi PLS ditulis sendiri (numpy/scipy) dan dibandingkan dengan paket `plspm`: selisih koefisien jalur ≤ 0,011 dan selisih loading ≤ 0,03.

## Etika dan anonimitas
Data tidak memuat identitas partisipan. Partisipasi bersifat sukarela dengan persetujuan sadar (lihat Pernyataan Persetujuan Etik pada artikel).

## Lisensi
Data/hasil: CC BY 4.0. Kode: MIT. Lihat `LICENSE`. Sitasi: lihat `CITATION.cff`.
