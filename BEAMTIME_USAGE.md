# Using AstraXAS Beamtime Mode at ASTRA

A step-by-step procedure for live QC during a real experiment.
Read the installation section before the beamtime; refer to the
operating sections during.

This guide covers Windows installation specifically, since that's
the platform on the ASTRA workstation. Linux notes are included
where the procedures differ.

---

## Part 1 — Installation on the ASTRA workstation (Windows)

Do this before the experiment starts. If possible, do it at least
a few hours in advance so you have time to debug any Windows-
specific issues without time pressure.

### What you need

- Windows 10 or Windows 11
- Administrator access OR ability to install Python for the
  current user (sufficient for our purposes — full admin not
  strictly required if Python and pip work)
- Internet access to download Python and pip packages
- About 2 GB of free disk space (Python + dependencies +
  a working area)
- The AstraXAS source code (either via `git clone` or a zip
  download from GitHub)

### Step 1: Install Python

If Python isn't already installed, get **Python 3.11 or 3.12** from
[python.org](https://www.python.org/downloads/windows/).

Choose the **64-bit installer** for your Windows version.

During the installer:

- ✅ **Check "Add Python to PATH"** at the bottom of the first
  screen. This is critical. If you forget this, you'll have to
  edit environment variables manually later.
- ✅ Choose "Install for all users" if you have admin rights.
  Otherwise "Install for current user only" works fine.
- Click "Install Now" with default settings.

After install, verify in a fresh PowerShell or Command Prompt
window:

```powershell
python --version
pip --version
```

You should see something like `Python 3.11.x` and `pip 24.x`. If
you get "command not found" errors, the PATH didn't get updated —
either re-run the installer with the PATH checkbox enabled, or
add `C:\Users\<you>\AppData\Local\Programs\Python\Python311\` and
its `Scripts\` subdirectory to PATH manually.

### Step 2: Install Git (if not already present)

The beamline workstation may already have Git. Check first:

```powershell
git --version
```

If you get a version number, skip ahead to Step 3.

If not, download Git for Windows from
[git-scm.com](https://git-scm.com/download/win). Use default
installer settings. After install, verify in a **new** PowerShell
window:

```powershell
git --version
```

### Step 3: Clone the AstraXAS repository

Pick a location to install AstraXAS. The user's home directory is
a good default. In PowerShell:

```powershell
cd $HOME
git clone https://github.com/cemcelikutku/AstraXAS.git
cd AstraXAS
```

You should now be inside the AstraXAS folder. Verify by listing
the contents:

```powershell
ls
```

You should see folders like `astra_xas`, `tests`, `examples`,
plus files like `README.md` and `requirements.txt`.

### Step 4: Create a virtual environment

This is **strongly recommended** so the AstraXAS dependencies
don't conflict with anything else installed on the workstation:

```powershell
python -m venv venv
```

Activate it (every new terminal session needs this):

```powershell
.\venv\Scripts\Activate.ps1
```

If you get an error about execution policy, run this once (it
allows scripts to run for your current user only):

```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```

Then activate again. Once activated, you should see `(venv)` at
the start of your prompt.

### Step 5: Install dependencies

With the venv activated:

```powershell
pip install --upgrade pip
pip install -r requirements.txt
pip install -r requirements-dev.txt
```

This installs `numpy`, `scipy`, `matplotlib`, `xraylarch`,
`watchdog`, `pyyaml`, and `pytest`. The `xraylarch` install in
particular can take a few minutes — it's a large package with
many dependencies.

If any install fails, the most common Windows-specific cause is
missing C++ build tools for some scientific package. The error
message will usually say "Microsoft Visual C++ 14.0 or greater is
required." If that happens, install the **Microsoft C++ Build
Tools** from
[visualstudio.microsoft.com/visual-cpp-build-tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/),
then retry the pip install.

### Step 6: Verify the install

Still in the activated venv:

```powershell
pytest
```

All tests should pass. You should see something like
`6 passed in 30.45s` (or similar — the count may have grown).

If any test fails, **fix it before the beamtime**, not during.
The most likely Windows-specific failure modes:

- **Filesystem watcher tests timing out**: Windows's filesystem
  events behave slightly differently from Linux. The Phase 1-3
  tests should already account for this, but if a test times
  out, it might be a race condition that needs the timeout
  raised on Windows.
- **Path separator issues**: The codebase uses `pathlib.Path`
  throughout, so this should be fine, but if you see errors
  involving backslashes vs forward slashes, file an issue.
- **PowerShell vs bash differences**: Some test fixtures use
  shell commands. They should all use Python primitives, but
  if you see bash-specific failures, file an issue.

### Step 7: Verify the CLI works

```powershell
python -m astra_xas.beamtime --help
python -m astra_xas.beamtime watch --help
python -m astra_xas.beamtime replay --help
```

All three should print help text without errors.

### Step 8: Do a dry run with synthetic data

Before relying on this at a real beamtime, do one full end-to-end
test with synthetic data. This catches Windows-specific issues
that the unit tests might miss.

Open **two PowerShell windows** (you'll need both later).

**In Terminal 1**, with the venv activated:

```powershell
cd $HOME\AstraXAS
.\venv\Scripts\Activate.ps1
python -m astra_xas.beamtime watch C:\temp\astra_smoke_test
```

This should start the watcher. It'll create
`C:\temp\astra_smoke_test\` and `C:\temp\astra_smoke_test-beamtime\`
automatically. Leave this terminal alone.

**In Terminal 2** (also activate the venv):

```powershell
cd $HOME\AstraXAS
.\venv\Scripts\Activate.ps1
```

Then drop two synthetic scans:

```powershell
python -c @"
from pathlib import Path
from astra_xas.beamtime._synthetic import write_synthetic_xasd
import shutil, time

watch_folder = Path(r'C:\temp\astra_smoke_test')
watch_folder.mkdir(parents=True, exist_ok=True)

for i in range(1, 3):
    src = Path(rf'C:\temp\smoke_{i}_source.xasd')
    final_dst = watch_folder / f'smoke_test_{i}.xasd'
    tmp_dst = watch_folder / f'smoke_test_{i}.xasd.tmp'
    write_synthetic_xasd(src, seed=i)
    shutil.copy2(src, tmp_dst)
    tmp_dst.replace(final_dst)
    print(f'wrote smoke_test_{i}.xasd')
    time.sleep(2)
"@
```

You should see two `wrote ...` lines in Terminal 2, and Terminal 1
should print two corresponding status lines.

Then open the dashboard. In Terminal 2:

```powershell
start C:\temp\astra_smoke_test-beamtime\index.html
```

This opens your default browser pointing at the dashboard. You
should see two scans under "Recent scans" and one group under
"Live groups" with `replicates: 2`, `status: ready`, and a
thumbnail.

If all of that works, **the system is fully functional on this
machine** and you can use it for the real beamtime.

In Terminal 1, press `Ctrl-C` to stop the watcher.

### Step 9: Prepare configs and log directory

```powershell
mkdir $HOME\astra_configs
mkdir $HOME\beamtime_logs
```

Then create at least one config JSON for the element you'll be
measuring. For example, for P K-edge, save the following as
`C:\Users\<you>\astra_configs\p_k.json`:

```json
{
  "analysis_mode": "fluo",
  "e0": 2145.5,
  "pre1": -60.0,
  "pre2": -15.0,
  "norm1": 40.0,
  "norm2": 250.0,
  "align_window_min": 2140.0,
  "align_window_max": 2160.0,
  "plot_energy_min": 2120.0,
  "plot_energy_max": 2300.0,
  "shift_bound_min": -3.0,
  "shift_bound_max": 3.0
}
```

You can edit this in Notepad, VS Code, or any text editor. Save
as plain text with `.json` extension. Do not save as `.json.txt`
— Windows sometimes hides the `.txt` and you'll spend 20 minutes
wondering why the file isn't being read.

If you'll measure multiple elements, create one config per
element (`s_k.json`, `cu_k.json`, etc.).

### Common Windows install pitfalls

- **"python is not recognized"** after install: PATH didn't update.
  Open a **new** terminal, or log out and back in.
- **`pip install` fails with SSL errors**: Beamline networks
  sometimes have weird proxies. Try
  `pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org -r requirements.txt`
  or ask the beamline staff about proxy settings.
- **Anti-virus blocking watchdog**: Some Windows endpoints flag
  filesystem watchers as suspicious. If the watcher seems to start
  but never sees new files, check Windows Defender's quarantine.
- **PowerShell scripts blocked**: Run
  `Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser`
  once.
- **OneDrive interfering**: If your Documents folder is synced to
  OneDrive, file events can have weird latency. Use a path outside
  OneDrive for the watch folder (e.g., `C:\beamtime\` rather than
  `C:\Users\<you>\Documents\beamtime\`).

---

## Part 2 — What this tool does and doesn't do

**It does:** Watch a folder, detect when new `.xasd` files arrive,
run per-scan QC on each, group replicates by base name, and
merge-normalize groups when ≥2 replicates exist. All of this
happens automatically in the background. You see the results in a
web browser tab that auto-refreshes every 5 seconds.

**It doesn't:** Replace the beamline acquisition software. Your
existing workflow for setting up scans, energy calibration, sample
positioning — none of that changes. AstraXAS just reads the files
your acquisition system writes.

**You interact with the tool through:**
1. A terminal where the watcher is running (mostly: leave it alone)
2. A browser tab with the dashboard (mostly: glance at it)
3. The output directory (occasionally: when something looks wrong)

You do not type commands at the watcher while it's running. There's
no interactive prompt. It's a "set it up at session start, glance at
the dashboard during, Ctrl-C at session end" tool.

---

## Part 3 — At the beamline — first 10 minutes

### Step 1: Activate the environment

Every new terminal session needs the venv activated. In PowerShell:

```powershell
cd $HOME\AstraXAS
.\venv\Scripts\Activate.ps1
```

You should see `(venv)` at the start of your prompt.

### Step 2: Find the beamline's data directory

The most important thing is finding **where the acquisition
software writes `.xasd` files**. Possibilities:

- A fixed local path: `C:\beamline_data\` or similar
- A network share mapped to a drive letter
- Your own user directory

Ask the beamline staff if you don't see it immediately. Once you
know the path, write it down. Let's call it `<DATA_DIR>` for the
rest of this guide.

**Important:** You want the directory where files are written *as
they're acquired*, not a directory you have to manually copy into.
The watcher needs to see new files appear in real time.

### Step 3: Verify the file naming convention

In the data directory, look at filenames from a previous experiment
(or run one quick test scan). The watcher groups replicates by
splitting the trailing `_<number>` from the filename:

- `mySample_1.xasd` → group `mySample`, replicate 1 ✓
- `mySample_2.xasd` → group `mySample`, replicate 2 ✓
- `mySample_a.xasd` → group `mySample_a`, no replicate ✗
- `mySample.xasd` → group `mySample`, no replicate ✗

If the beamline software writes filenames in a format that doesn't
end in `_<digit>`, you have two options:
- (a) Configure the acquisition software to use the convention
  if possible.
- (b) Live with each scan being its own group (you still get
  per-scan QC, just no group merge-normalize).

### Step 4: Open two PowerShell windows

You need two terminals:
- **Terminal 1:** the watcher will run here. Set it up once,
  leave it alone until session end.
- **Terminal 2:** for everything else.

Both terminals need the venv activated.

### Step 5: Start the watcher

In **Terminal 1**:

```powershell
python -m astra_xas.beamtime watch <DATA_DIR> `
    -c $HOME\astra_configs\p_k.json `
    -l $HOME\beamtime_logs\$(Get-Date -Format 'yyyy-MM-dd')_session1.log
```

(The backticks in PowerShell are line-continuation characters,
equivalent to `\` in bash. You can also put it all on one line.)

Replace `<DATA_DIR>` with the actual path, and `p_k.json` with
the config matching your element.

You should see a single line:

```
[2026-05-15T09:14:33] AstraXAS Beamtime watcher starting; tee log -> C:\Users\you\beamtime_logs\2026-05-15_session1.log
```

The watcher is now running. **Do not close this terminal.** Do
not type into it. Leave it visible if you can.

### Step 6: Open the dashboard

In **Terminal 2**:

```powershell
start <DATA_DIR>-beamtime\index.html
```

This opens your default browser pointing at the dashboard.
**Leave this tab open.** It auto-refreshes every 5 seconds.

### Step 7: Verify with one test scan

Before doing science, sanity-check end-to-end with one real scan.
Have the beamline software acquire a single `.xasd` file. Watch
what happens in three places:

1. **Terminal 1:** within ~2 seconds the watcher should print a
   status line.
2. **Browser dashboard:** within 5 seconds (next auto-refresh),
   you should see one row appear under "Recent scans".
3. **Log file:** if you `Get-Content` the log file in Terminal 2
   (`Get-Content $HOME\beamtime_logs\<file>.log`), you should
   see the same line.

If all three show the scan, you're good to start the experiment.
**If any of them doesn't show it, stop and debug now** — see the
troubleshooting section below.

---

## Part 4 — During the experiment — the rhythm

Once everything is running, you mostly leave the watcher alone.

### What you do actively

- **Run your experiment normally.** AstraXAS reads files; it
  doesn't acquire them.
- **Glance at the dashboard tab between scan groups.** Look at:
  - Are scans appearing as expected?
  - What status does each scan have?
  - Once you have ≥2 replicates of the same group, the "Live
    groups" section should appear with `ready` status and a
    thumbnail of the merged spectrum.
- **Look at the group QC plot when a group goes ready.** Click
  the thumbnail to see full size. You're looking for: do the
  replicates overlay cleanly? Is the merged spectrum
  recognizable as the element you're measuring?

### What you do NOT do

- Don't restart the watcher between samples. It handles multiple
  groups simultaneously.
- Don't manually edit files in the output directory while the
  watcher is running.
- Don't switch configs mid-session unless you're switching
  elements (see below).

### Switching elements mid-session

If you change to a different element (e.g., Cu K after P K), the
existing watcher's config no longer applies.

**Recommended approach:** Stop the watcher (`Ctrl-C` in
Terminal 1), restart with the new config:

```powershell
python -m astra_xas.beamtime watch <DATA_DIR> `
    -c $HOME\astra_configs\cu_k.json `
    -l $HOME\beamtime_logs\$(Get-Date -Format 'yyyy-MM-dd')_session2.log
```

The previous element's groups stay in the dashboard (they're
persisted to disk), and new scans use the new config. Use a
different log file name so the two sessions don't get mixed up.

### When something looks wrong

**Lots of `warn` status with high jump counts.** Usually fine.
The detector-jump heuristic flags real instrumental noise that
doesn't affect the science. If the merged spectrum looks right,
ignore the warn.

**A scan rejected.** Check the session log:

```powershell
Get-Content <DATA_DIR>-beamtime\ASTRA_beamtime_session.log
```

The `notes` column tells you why. Common reasons:
- `load_failed:` — file couldn't be parsed. Could be a partial
  write the watcher caught too early, an actual file format
  problem, or wrong file extension.
- `pipeline_failed:` — `_entry_from_scan` failed. Usually a
  column mismatch or unexpected file structure.

**A group's status is `error`.** Check the JSON summary:

```powershell
Get-Content <DATA_DIR>-beamtime\groups\<group_name>_group_summary.json
```

The `last_merge_error` field has the exception. Most likely
cause: the config's energy windows don't fit inside your data
range.

**Dashboard not refreshing.** Browser tab may be throttled.
Click into it or hit F5.

---

## Part 5 — End of session — clean shutdown

### Step 1: Stop the watcher

In Terminal 1, press `Ctrl-C`. The watcher should exit within
~5 seconds.

### Step 2: Verify outputs

```powershell
ls <DATA_DIR>-beamtime\groups\
ls <DATA_DIR>-beamtime\plots\group_qc\
(Get-Content <DATA_DIR>-beamtime\ASTRA_beamtime_session.log).Count
```

Expected:
- One `_processed.dat`, `_norm.dat`, `_flat.dat`, and
  `_group_summary.json` per group
- One `_replicate_qc.png` per group
- The session log row count should be roughly the number of
  scans you took plus a few overhead rows

### Step 3: Backup everything

This is the most important step. Beamline workstations get reset
between users.

To a USB drive (say it's mounted as `E:\`):

```powershell
Copy-Item -Recurse <DATA_DIR>-beamtime E:\beamtime_$(Get-Date -Format 'yyyy-MM-dd')\
Copy-Item $HOME\beamtime_logs\*.log E:\beamtime_$(Get-Date -Format 'yyyy-MM-dd')\
```

Or to your own laptop via a network share, scp, or whatever's
available.

### Step 4: Take notes for after

The single most valuable thing from this beamtime is **the list
of things AstraXAS got wrong, missed, or made awkward**. Write
these down while they're fresh:

- What worked smoothly?
- What surprised you?
- What did you wish was there but wasn't?
- What was hard that should be easy?
- What looked wrong but turned out to be fine, or vice versa?

---

## Part 6 — Troubleshooting reference

### "Watcher prints nothing when scans arrive"

The watcher might be watching the wrong directory, or the
acquisition software is writing to a different path than you
think.

```powershell
ls <DATA_DIR>\*.xasd
```

If this lists files but the watcher hasn't reacted, double-check
the path matches exactly. PowerShell is case-insensitive but
some Windows filesystems aren't entirely consistent — use
copy-paste rather than retyping.

### "Watcher starts but immediately exits"

Probably a JSON syntax error in the config file. Verify:

```powershell
python -c "import json; json.load(open(r'C:\Users\you\astra_configs\p_k.json'))"
```

If that fails, fix the JSON. Common Windows mistakes: trailing
commas, single quotes instead of double quotes, BOM at start of
file (saved with the wrong encoding in Notepad).

### "All scans get rejected"

Run one scan through the loader to see if the data is parseable:

```powershell
python -c @"
from astra_xas.io import load_xasd
from pathlib import Path
data = load_xasd(Path(r'<DATA_DIR>\some_scan.xasd'))
print(list(data.keys()))
"@
```

If this raises an error, the file format isn't what `load_xasd`
expects.

### "Group QC plot looks wrong"

Open the merged `.dat` files:

```powershell
Get-Content <DATA_DIR>-beamtime\groups\<group>_norm.dat -Head 20
```

If the `norm` column is all near-zero or all-NaN, the pre-edge
or post-edge windows in your config don't fit your data. Stop
the watcher, edit the config, restart. The next replicate will
trigger a re-merge.

### "Dashboard shows scans but no Live groups section"

You need ≥2 accepted replicates of the same `base_name`. Check:

```powershell
Get-Content <DATA_DIR>-beamtime\ASTRA_beamtime_session.log
```

Are filenames following `<base>_<N>.xasd`? Are at least 2
accepted (not reject)?

### "WindowsError or PermissionError on file operations"

OneDrive, Dropbox, or anti-virus may be holding files open. Try
moving the watch folder to a path that isn't synced or
monitored.

### "Watchdog stops detecting changes after a while"

Rare but possible on Windows with long-running watchers. If you
notice the dashboard stops updating despite scans appearing in
the data folder, restart the watcher. The group state persists
across restarts.

---

## What this guide assumes

- Working AstraXAS install on the ASTRA workstation (Part 1
  covers this)
- At least one config JSON for the element you're measuring
- You can identify where the beamline software writes files
- You have permission to read that directory and write to a
  chosen output location

If any of these isn't true, fix it before the experiment starts.
The beamtime itself is not the time to debug installation, paths,
or permissions.
