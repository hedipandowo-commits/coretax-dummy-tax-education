# Kamus data: `data_respons_individual.csv`

96 baris (satu baris = satu partisipan), 43 kolom. Tidak ada nilai kosong. Data **anonim**: `ID` hanya nomor urut (R001–R096); tidak ada nama, NIM, e-mail, institusi, atau data pribadi lain.

| Kolom | Jumlah | Keterangan | Skala |
|---|---|---|---|
| `ID` | 1 | Kode partisipan anonim | R001–R096 |
| `PU1`–`PU6` | 6 | Perceived Usefulness (persepsi kegunaan) | Likert 1–5 |
| `PEOU1`–`PEOU6` | 6 | Perceived Ease of Use (persepsi kemudahan) | Likert 1–5 |
| `ATU1`–`ATU6` | 6 | Attitude Toward Using (sikap) | Likert 1–5 |
| `BI1`–`BI6` | 6 | Behavioral Intention (niat berperilaku) | Likert 1–5 |
| `TCSE1`–`TCSE6` | 6 | Tax Compliance Self-Efficacy (butir pada Lampiran 1 artikel) | Likert 1–5 |
| `SUS1`–`SUS10` | 10 | System Usability Scale (butir ganjil positif, genap negatif) | Likert 1–5 |
| `PRE` | 1 | Skor mentah pre-test | 0–30 (30 soal pilihan ganda) |
| `POST` | 1 | Skor mentah post-test | 0–30 |

## Pengolahan (dilakukan oleh skrip, data mentah tidak diubah)
- Skor pre/post dikonversi ke skala 0–100: `skor / 30 × 100`.
- Skor SUS: butir ganjil `(x−1)`, butir genap `(5−x)`, dijumlahkan lalu `× 2,5` (rentang 0–100).
- Skor konstruk PU, PEOU, ATU, BI, TCSE = rerata enam butir.
- N-Gain = (post − pre) / (100 − pre).
