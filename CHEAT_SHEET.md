> **READ BEAMTIME_USAGE.md ONCE BEFORE THE EXPERIMENT.** This sheet
> is a reference, not a tutorial.

# AstraXAS beamtime cheat sheet

Quick reference for AstraXAS at the beamtime. See
[BEAMTIME_USAGE.md](BEAMTIME_USAGE.md) for full installation and
operating procedure.

Assumes Windows + PowerShell + venv already set up at
`C:\Users\cemcl\AstraXAS`. Replace example paths below with the
actual beamline paths on the day.

Example paths used throughout:

- `DATA_DIR    = C:\Users\cemcl\beamtime_data`
- `CONFIG      = C:\Users\cemcl\AstraXAS\configs\k_k_operando_pd_foil.json`
- `LOG         = C:\Users\cemcl\beamtime_logs\session1.log`
- `OUTPUT_FINAL = C:\Users\cemcl\beamtime_data-final`

## 1. The 5 commands

### 1. Activate venv (every new terminal)

```powershell
cd $HOME\AstraXAS
.\venv\Scripts\Activate.ps1
```

Verify the prompt now starts with `(venv)`. If not, stop and fix
before doing anything else.

### 2. Start the watcher — Terminal 1

```powershell
astra-xas-beamtime watch C:\Users\cemcl\beamtime_data `
    -c C:\Users\cemcl\AstraXAS\configs\k_k_operando_pd_foil.json `
    -l C:\Users\cemcl\beamtime_logs\session1.log
```

Leave this terminal alone. Do not type into it.

### 3. Open the dashboard — Terminal 2

```powershell
cd $HOME\AstraXAS
.\venv\Scripts\Activate.ps1
start C:\Users\cemcl\beamtime_data-beamtime\index.html
```

Browser tab auto-refreshes every 5 seconds. Leave it open.

### 4. Stop the watcher

`Ctrl-C` in Terminal 1. Expected output:

```
Watcher stopped cleanly. Total scans processed: N
```

### 5. Run the offline pipeline — AFTER the experiment

```powershell
astra-xas C:\Users\cemcl\beamtime_data `
    --config C:\Users\cemcl\AstraXAS\configs\k_k_operando_pd_foil.json `
    -o C:\Users\cemcl\beamtime_data-final
```

**This is the run that produces your scientific results. Live mode
is QC only.**

## 2. What live mode does and doesn't do

**DOES (during the experiment):**

- Validates each scan, computes processed μ(E), runs detector jump
  diagnostics
- Writes per-scan QC plots and the live HTML dashboard
- When 2+ replicates of the same group arrive: live merge + normalize

**DOES NOT (deferred to offline run):**

- Foil drift correction (no Pd foil shifts applied live)
- Cross-scan alignment between replicates
- Full validation / self-absorption diagnostic
- Comprehensive outlier rejection

**One-line takeaway: Live dashboard = QC. Offline pipeline = results.**

## 3. What to look for during the experiment

**GOOD signs:**

- Sample scans show K K-edge rise around 3608 eV in per-scan plot
- Pd foil scans appear every ~3 sample scans, show Pd L₁ at ~3604 eV
- Detector jump warnings stay rare
- Live groups dashboard fills in as replicates accumulate

**WARNING signs (investigate):**

- Many scans rejected → check validation reasons in session log
- No Pd foil scans appearing → verify `foil_keyword` matches actual
  filenames; if not, stop watcher, edit config, restart
- Detector warnings on every scan → beamline instability, ask staff

**EXPECTED but not problems:**

- Live group plots may look slightly mismatched between replicates
  (drift correction is offline)
- Energy axis may show edge slightly shifted (drift will be
  corrected offline)
- DO NOT change anything mid-experiment based on live merge appearance

## 4. If something goes wrong

- Session log: `<DATA_DIR>-beamtime\ASTRA_beamtime_session.log` —
  full per-scan trace
- Per-scan plots: `<DATA_DIR>-beamtime\plots\beamtime\<scan>.png`
- Tee log (if `-l` was used): human-readable transcript
- Worst case: `Ctrl-C`, fix config, restart watcher. Data already
  written is preserved.

## 5. End of experiment workflow

1. Stop watcher (`Ctrl-C` in Terminal 1)
2. Backup raw data (USB drive or network)
3. Run offline pipeline (section 1, command 5)

## 6. Pre-experiment verification

Before scans start arriving:

- [ ] venv activated (prompt shows `(venv)`)
- [ ] `astra-xas-beamtime --help` works
- [ ] Watch folder exists and is empty (or you've planned how to
      handle pre-existing files)
- [ ] Dashboard URL opens correctly with synthetic test scan
