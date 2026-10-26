# ViHoP demo (Gossenkollesee selected cells)

Real per-cell combined FASTQs from the Gossenkollesee analysis (two barcodes), plus the
host annotations and viral lastal hits needed to run the host–virus pairing
summary stage on a laptop. Sample IDs below are the original 10x library names.

| Cell | Sample ID | Combined FASTQ |
|------|-----------|----------------|
| `AATCGGCTCTGGATGC` | Alps5 | ~850 KB |
| `ACATCCCTCTATGCCC` | Alps_Amir | ~620 KB |

Also included: `transcripts.edit.fasta` for each cell (rnaSPAdes products from the run).

## Run

```bash
# from the ViHoP root, with the vihop conda env active
bash demo/run_demo.sh
```

**Expected runtime:** under 1 minute on a normal desktop/laptop.  
**Expected result:** script exits 0 and prints `PASS`; outputs match `demo/expected/`.

This demo does **not** re-run Cell Ranger or assembly. Full Gossenkollesee reproduction uses `config/pairing.config.Gossenkollesee` on an LSF cluster (see main README).
