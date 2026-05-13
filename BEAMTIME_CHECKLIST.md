# Beamtime Checklist

A practical procedure for using AstraXAS Beamtime Mode during a
synchrotron experiment. Designed to be readable when tired.

## Day before the beamtime

- [ ] `git pull` latest AstraXAS on the laptop you're bringing
- [ ] `pytest` — all phases should pass without errors
- [ ] `python -m astra_xas.beamtime --help` works without errors
- [ ] Confirm `xraylarch` and `watchdog` are installed in the
      Python environment you'll use
- [ ] Identify which element(s) you'll be measuring and prepare
      a config JSON for each (see `configs/`)
- [ ] Verify each config has correct `e0`, `align_window_min`,
      `align_window_max`, `pre1`/`pre2`/`norm1`/`norm2`, and
      `analysis_mode`
- [ ] Charge laptop battery, pack charger
- [ ] Bring a USB drive for backing up the output directory

## Arriving at the beamline

- [ ] Find the directory where the beamline software writes
      `.xasd` files. Confirm the filename pattern is
      `<sample_name>_<replicate_number>.xasd` (so
      `split_replicate_suffix` parses it correctly)
- [ ] Choose where AstraXAS output should go. By default the
      watcher writes to `<incoming_dir>-beamtime/`, which is
      usually fine
- [ ] Open one terminal for the watcher, one for ad-hoc
      commands. Don't multiplex everything into one window
- [ ] Have the matching config JSON ready in a known path

## Starting a session

```bash
cd ~/Applications/astra-xas-processor/astra-xas-processor

python -m astra_xas.beamtime watch /path/to/beamline/incoming \
    -c ~/configs/<element>.json \
    -l ~/beamtime_logs/$(date +%F).log
```

- [ ] Watcher prints "starting" line and stays running
- [ ] Open the dashboard in a browser tab and leave it open:
      `file:///path/to/beamline/incoming-beamtime/index.html`
- [ ] Verify the dashboard loads with empty Recent scans /
      no Live groups (correct empty state)
- [ ] Drop one test scan into the incoming folder. Verify it
      appears in the watcher's stdout, on the dashboard, and
      in the tee log file. If any of these three doesn't show
      up, stop and debug before doing real measurements.

## During the session

- [ ] Glance at the dashboard between scan groups. Specifically
      look at:
      - Group QC plots: are replicates overlaying cleanly?
      - Per-scan plots: do the raw detector channels look sane
        (no clipping, no all-zero channels)?
      - Status colors: any reject? Why?
- [ ] If something looks wrong, check the tee log file
      (`tail -f` is your friend)
- [ ] If you change samples or elements, decide whether to:
      (a) Stop the watcher and restart with a different config
          (groups from previous element get persisted as-is and
          will reappear on dashboard if you point at the same
          output_dir)
      (b) Keep the same watcher running if the same config
          parameters apply (rare)

## End of session

- [ ] In the watcher's terminal: `Ctrl-C`. Watcher should exit
      within ~5 seconds
- [ ] Verify session log has the expected number of rows:
      `wc -l <output_dir>/ASTRA_beamtime_session.log`
- [ ] Backup the entire output directory to USB:
      `cp -r <output_dir> /media/usb/<beamtime_id>-<date>/`
- [ ] Backup the tee log file too: `cp ~/beamtime_logs/*.log
      /media/usb/<beamtime_id>-<date>/`

## After the beamtime — feedback capture

The most valuable output of this trip isn't the data, it's the
list of things AstraXAS got wrong, missing features you wished
existed, and UX papercuts.

- [ ] Open a new file `feedback/<date>.md` in the repo and
      write down everything that was annoying, unexpected, or
      missing
- [ ] Note which detector-jump thresholds caused too many
      warns (e.g., "P K at 10.0 was too sensitive — try 25")
- [ ] Note any file-naming patterns the beamline uses that
      don't match `split_replicate_suffix`
- [ ] Note dashboard features you wished you had (live foil
      drift indicator? sample-name search? per-element preset
      switcher?)

## Common problems and fixes

### "Watcher doesn't see my scans"
- Verify the path passed to `watch` is the actual directory
  where the beamline writes files
- Verify the file extension is exactly `.xasd` (not `.XASD` or
  `.xasd.bak`)
- Check the watcher's stdout for "WARNING: file did not settle"
  messages (means the file wasn't fully written before the
  watcher tried to process it)

### "Status is reject for everything"
- Either `load_xasd` is failing (file format mismatch) or
  `_validate_processing_inputs` is failing (config doesn't
  match the data, e.g., wrong analysis_mode, energy range
  outside data bounds)
- Check tee log for the `notes` field on a reject row to see
  the specific error

### "Group never becomes ready"
- Group needs >= 2 accepted (ok or warn) replicates with the
  same `base_name`
- Check that filenames follow `<base>_<N>.xasd` pattern (e.g.,
  `mySample_1.xasd`, `mySample_2.xasd`)
- Check session log: were both replicates accepted (status ok
  or warn) or was one rejected?

### "Dashboard shows broken image icons"
- Probably a path encoding issue with unusual characters in
  filenames. Open `index.html` in the browser's developer
  tools and check the `img src` URLs
- File a `feedback/` note for this — we'd want to fix it
