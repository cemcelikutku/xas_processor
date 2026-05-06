from __future__ import annotations

import json
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

import numpy as np
from scipy.signal import savgol_filter


def _safe_savgol(y: np.ndarray, window_length: int, polyorder: int) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 3:
        return y.copy()
    if window_length % 2 == 0:
        window_length += 1
    if window_length >= n:
        window_length = n - 1 if n % 2 == 0 else n
    if window_length % 2 == 0:
        window_length -= 1
    if window_length < 3:
        return y.copy()
    if polyorder >= window_length:
        polyorder = window_length - 1
    if polyorder < 1:
        return y.copy()
    return savgol_filter(y, window_length=window_length, polyorder=polyorder)


def read_dat_table(path: str | Path) -> tuple[np.ndarray, dict[str, np.ndarray], list[str]]:
    """Read an ASTRA .dat file with two or more numeric columns.

    Handles both processed outputs (energy + one signal) and detector_raw
    exports (energy + I0/I1/I2/IF/FDT/Ir + derived signals).
    """
    path = Path(path)
    if not path.exists():
        raise ValueError(f"File does not exist: {path}")

    header: list[str] | None = None
    rows: list[list[float]] = []

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line_no, line in enumerate(fh, start=1):
            text = line.strip()
            if not text:
                continue

            # Try to recover column names from commented headers too.
            if text.startswith("#"):
                clean = text.lstrip("#").strip()
                if clean.lower().startswith("columns:"):
                    header = clean.split(":", 1)[1].strip().split()
                elif clean.split() and clean.split()[0].lower() in {"energy_ev", "energy", "e", "energy/ev"}:
                    header = clean.split()
                continue

            parts = text.split()
            if len(parts) < 2:
                continue

            try:
                values = [float(x) for x in parts]
            except ValueError:
                # Plain header line, e.g. "energy_eV I0 I1 I2 IF".
                header = parts
                continue

            rows.append(values)

    if not rows:
        raise ValueError(f"No numeric data rows found in {path.name}")

    ncols = max(set(len(r) for r in rows), key=[len(r) for r in rows].count)
    rows = [r for r in rows if len(r) == ncols]
    arr = np.asarray(rows, dtype=float)
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise ValueError(f"Expected at least two numeric columns in {path.name}")

    if not header or len(header) != arr.shape[1]:
        # Fallback for older two-column processed files.
        stem = path.stem.lower()
        if arr.shape[1] == 2:
            if stem.endswith("_norm"):
                header = ["energy_eV", "norm"]
            elif stem.endswith("_raw"):
                header = ["energy_eV", "raw"]
            elif stem.endswith("_flat"):
                header = ["energy_eV", "flat"]
            else:
                header = ["energy_eV", "signal"]
        else:
            header = ["energy_eV"] + [f"col{i}" for i in range(1, arr.shape[1])]

    lowered = [h.lower() for h in header]
    energy_idx = 0
    for i, name in enumerate(lowered):
        if name in {"energy_ev", "energy", "e", "energy/ev"}:
            energy_idx = i
            break

    energy = np.asarray(arr[:, energy_idx], dtype=float)
    columns: dict[str, np.ndarray] = {}
    for i, name in enumerate(header):
        if i == energy_idx:
            continue
        columns[name] = np.asarray(arr[:, i], dtype=float)

    mask = np.isfinite(energy)
    for values in columns.values():
        mask &= np.isfinite(values)
    if mask.sum() == 0:
        raise ValueError(f"No finite numeric data found in {path.name}")

    order = np.argsort(energy[mask])
    energy = energy[mask][order]
    clean_columns = {name: values[mask][order] for name, values in columns.items()}
    return energy, clean_columns, header


def _find_column(columns: dict[str, np.ndarray], *names: str) -> np.ndarray | None:
    exact = {k.lower(): k for k in columns}
    for name in names:
        key = name.lower()
        if key in exact:
            return columns[exact[key]]
    return None


def read_spectrum_dat(path: str | Path, channel: str = "auto") -> tuple[np.ndarray, np.ndarray, str]:
    energy, columns, _header = read_dat_table(path)
    if not columns:
        raise ValueError(f"No signal columns found in {Path(path).name}")

    channel = (channel or "auto").strip()
    channel_lower = channel.lower()

    if channel_lower in {"auto", "first signal", "first_signal"}:
        name = next(iter(columns))
        return energy, columns[name], name

    # Exact column names first.
    col = _find_column(columns, channel)
    if col is not None:
        return energy, col, channel

    # Friendly aliases and computed signals for detector raw files.
    if channel_lower in {"i0", "i1", "i2", "if", "fdt", "ir"}:
        col = _find_column(columns, channel)
        if col is None:
            raise ValueError(f"Column {channel!r} not found in {Path(path).name}")
        return energy, col, channel

    i0 = _find_column(columns, "I0")
    i1 = _find_column(columns, "I1")
    i2 = _find_column(columns, "I2")
    iff = _find_column(columns, "IF")

    eps = 1e-30
    if channel_lower in {"if/i0", "fluo", "mu_fluo", "mu_fluo_ifi0"}:
        existing = _find_column(columns, "mu_fluo_IFI0")
        if existing is not None:
            return energy, existing, "mu_fluo_IFI0"
        if iff is None or i0 is None:
            raise ValueError(f"IF/I0 requires IF and I0 columns in {Path(path).name}")
        return energy, iff / np.clip(i0, eps, None), "IF/I0"

    if channel_lower in {"ln(i0/i1)", "trans", "mu_trans", "mu_trans_lni0i1"}:
        existing = _find_column(columns, "mu_trans_lnI0I1")
        if existing is not None:
            return energy, existing, "mu_trans_lnI0I1"
        if i0 is None or i1 is None:
            raise ValueError(f"ln(I0/I1) requires I0 and I1 columns in {Path(path).name}")
        return energy, np.log(np.clip(i0, eps, None) / np.clip(i1, eps, None)), "ln(I0/I1)"

    if channel_lower in {"ln(i1/i2)", "ref", "mu_ref", "mu_ref_lni1i2"}:
        existing = _find_column(columns, "mu_ref_lnI1I2")
        if existing is not None:
            return energy, existing, "mu_ref_lnI1I2"
        if i1 is None or i2 is None:
            raise ValueError(f"ln(I1/I2) requires I1 and I2 columns in {Path(path).name}")
        return energy, np.log(np.clip(i1, eps, None) / np.clip(i2, eps, None)), "ln(I1/I2)"

    available = ", ".join(columns.keys())
    raise ValueError(f"Channel {channel!r} not available in {Path(path).name}. Available: {available}")


def read_two_column_dat(path: str | Path) -> tuple[np.ndarray, np.ndarray, str]:
    """Backward-compatible reader used by older viewer code."""
    return read_spectrum_dat(path, channel="auto")


class SpectrumViewer(tk.Toplevel):
    """Interactive viewer for processed ASTRA .dat spectra."""

    def __init__(self, master=None):
        super().__init__(master)
        self.title("ASTRA Spectrum Viewer")
        self.geometry("1180x760")
        self.minsize(980, 650)

        self.files: list[Path] = []
        self.labels: list[str] = []

        self.e_min = tk.StringVar(value="7100.0")
        self.e_max = tk.StringVar(value="7160.0")
        self.use_smoothing = tk.BooleanVar(value=True)
        self.show_raw_too = tk.BooleanVar(value=False)
        self.plot_channel = tk.StringVar(value="auto")
        self.sg_window = tk.StringVar(value="11")
        self.sg_polyorder = tk.StringVar(value="2")
        self.linewidth = tk.StringVar(value="3.0")
        self.legend_fontsize = tk.StringVar(value="12")
        self.fig_width = tk.StringVar(value="8.0")
        self.fig_height = tk.StringVar(value="5.5")
        self.xlabel = tk.StringVar(value="Energy (eV)")
        self.ylabel = tk.StringVar(value="Normalized XANES")
        self.legend_location = tk.StringVar(value="lower right")
        self.title_text = tk.StringVar(value="")
        self.grid_on = tk.BooleanVar(value=False)
        self.selected_label = tk.StringVar(value="")
        self.selected_file = tk.StringVar(value="No spectrum selected")
        self.status = tk.StringVar(value="Select processed .dat files to compare.")

        self._build()

    def _build(self):
        root = ttk.Frame(self, padding=10)
        root.pack(fill="both", expand=True)
        root.columnconfigure(0, weight=1)
        root.columnconfigure(1, weight=0)
        root.rowconfigure(1, weight=1)

        title = ttk.Label(root, text="Spectrum Viewer", font=("TkDefaultFont", 16, "bold"))
        title.grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 8))

        file_frame = ttk.LabelFrame(root, text="Spectra", padding=10)
        file_frame.grid(row=1, column=0, sticky="nsew", padx=(0, 10))
        file_frame.rowconfigure(0, weight=1)
        file_frame.columnconfigure(0, weight=1)

        columns = ("file", "label")
        self.tree = ttk.Treeview(file_frame, columns=columns, show="headings", selectmode="extended")
        self.tree.heading("file", text="File")
        self.tree.heading("label", text="Plot label")
        self.tree.column("file", width=390, stretch=True)
        self.tree.column("label", width=420, stretch=True)
        self.tree.grid(row=0, column=0, sticky="nsew")
        scroll = ttk.Scrollbar(file_frame, orient="vertical", command=self.tree.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        self.tree.configure(yscrollcommand=scroll.set)
        self.tree.bind("<<TreeviewSelect>>", self._on_select)
        self.tree.bind("<Double-1>", lambda _e: self._focus_label_entry())

        file_buttons = ttk.Frame(file_frame)
        file_buttons.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(10, 0))
        
        for i in range(4):
            file_buttons.columnconfigure(i, weight=1)

        ttk.Button(file_buttons, text="Add .dat files", command=self.add_files).grid(row=0, column=0, sticky="ew", padx=2, pady=2)
        ttk.Button(file_buttons, text="Remove", command=self.remove_selected).grid(row=0, column=1, sticky="ew", padx=2, pady=2)
        ttk.Button(file_buttons, text="Clear", command=self.clear_files).grid(row=0, column=2, sticky="ew", padx=2, pady=2)
        ttk.Button(file_buttons, text="Auto labels", command=self.auto_labels).grid(row=0, column=3, sticky="ew", padx=2, pady=2)

        ttk.Button(file_buttons, text="Move up", command=self.move_up).grid(row=1, column=0, sticky="ew", padx=2, pady=2)
        ttk.Button(file_buttons, text="Move down", command=self.move_down).grid(row=1, column=1, sticky="ew", padx=2, pady=2)
        ttk.Button(file_buttons, text="Sort A-Z", command=self.sort_by_filename).grid(row=1, column=2, columnspan=2, sticky="ew", padx=2, pady=2)

        selection_frame = ttk.Frame(file_buttons)
        selection_frame.grid(row=2, column=0, columnspan=4, sticky="ew", pady=(5, 0))
        
        for i in range(3):
            selection_frame.columnconfigure(i, weight=1)
            
        ttk.Button(selection_frame, text="Select all", command=self.select_all).grid(row=0, column=0, sticky="ew", padx=2)
        ttk.Button(selection_frame, text="Deselect all", command=self.deselect_all).grid(row=0, column=1, sticky="ew", padx=2)
        ttk.Button(selection_frame, text="Invert selection", command=self.invert_selection).grid(row=0, column=2, sticky="ew", padx=2)

        editor = ttk.LabelFrame(file_frame, text="Selected spectrum label", padding=10)
        editor.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(10, 0))
        editor.columnconfigure(1, weight=1)
        ttk.Label(editor, text="File").grid(row=0, column=0, sticky="w", padx=(0, 8))
        ttk.Label(editor, textvariable=self.selected_file, foreground="#444").grid(row=0, column=1, sticky="w")
        ttk.Label(editor, text="Plot label").grid(row=1, column=0, sticky="w", padx=(0, 8), pady=(6, 0))
        self.label_entry = ttk.Entry(editor, textvariable=self.selected_label)
        self.label_entry.grid(row=1, column=1, sticky="ew", pady=(6, 0))
        self.label_entry.bind("<Return>", lambda _e: self.apply_label())
        ttk.Button(editor, text="Apply label", command=self.apply_label).grid(row=1, column=2, sticky="ew", padx=(8, 0), pady=(6, 0))

        session_buttons = ttk.Frame(file_frame)
        session_buttons.grid(row=3, column=0, columnspan=2, sticky="ew", pady=(8, 0))
        session_buttons.columnconfigure(0, weight=1)
        session_buttons.columnconfigure(1, weight=1)
        ttk.Button(session_buttons, text="Save viewer session", command=self.save_session).grid(row=0, column=0, sticky="ew", padx=(0, 4))
        ttk.Button(session_buttons, text="Load viewer session", command=self.load_session).grid(row=0, column=1, sticky="ew", padx=(4, 0))

        settings = ttk.LabelFrame(root, text="Plot settings", padding=10)
        settings.grid(row=1, column=1, sticky="nsew")
        settings.columnconfigure(1, weight=1)

        r = 0
        ttk.Label(settings, text="Energy range / eV").grid(row=r, column=0, sticky="w")
        range_frame = ttk.Frame(settings)
        range_frame.grid(row=r, column=1, sticky="w", pady=3)
        ttk.Entry(range_frame, textvariable=self.e_min, width=9).grid(row=0, column=0)
        ttk.Label(range_frame, text="to").grid(row=0, column=1, padx=4)
        ttk.Entry(range_frame, textvariable=self.e_max, width=9).grid(row=0, column=2)
        r += 1

        ttk.Label(settings, text="Y channel").grid(row=r, column=0, sticky="w", pady=3)
        ttk.Combobox(
            settings,
            textvariable=self.plot_channel,
            width=18,
            state="readonly",
            values=(
                "auto", "I0", "I1", "I2", "IF", "FDT", "Ir",
                "IF/I0", "ln(I0/I1)", "ln(I1/I2)",
                "mu_fluo_IFI0", "mu_trans_lnI0I1", "mu_ref_lnI1I2",
            ),
        ).grid(row=r, column=1, sticky="w", pady=3)
        r += 1

        ttk.Checkbutton(settings, text="Savitzky-Golay smoothing", variable=self.use_smoothing).grid(row=r, column=0, columnspan=2, sticky="w", pady=3)
        r += 1
        ttk.Checkbutton(settings, text="Show raw too", variable=self.show_raw_too).grid(row=r, column=0, columnspan=2, sticky="w", pady=3)
        r += 1
        ttk.Checkbutton(settings, text="Show grid", variable=self.grid_on).grid(row=r, column=0, columnspan=2, sticky="w", pady=3)
        r += 1

        for label, var in [
            ("SG window", self.sg_window),
            ("SG polyorder", self.sg_polyorder),
            ("Line width", self.linewidth),
            ("Legend font size", self.legend_fontsize),
        ]:
            ttk.Label(settings, text=label).grid(row=r, column=0, sticky="w", pady=3)
            ttk.Entry(settings, textvariable=var, width=10).grid(row=r, column=1, sticky="w", pady=3)
            r += 1

        ttk.Label(settings, text="Figure size").grid(row=r, column=0, sticky="w", pady=3)
        size_frame = ttk.Frame(settings)
        size_frame.grid(row=r, column=1, sticky="w", pady=3)
        ttk.Entry(size_frame, textvariable=self.fig_width, width=7).grid(row=0, column=0)
        ttk.Label(size_frame, text="x").grid(row=0, column=1, padx=3)
        ttk.Entry(size_frame, textvariable=self.fig_height, width=7).grid(row=0, column=2)
        r += 1

        ttk.Label(settings, text="Plot title").grid(row=r, column=0, sticky="w", pady=3)
        ttk.Entry(settings, textvariable=self.title_text, width=24).grid(row=r, column=1, sticky="ew", pady=3)
        r += 1
        ttk.Label(settings, text="X-axis label").grid(row=r, column=0, sticky="w", pady=3)
        ttk.Entry(settings, textvariable=self.xlabel, width=24).grid(row=r, column=1, sticky="ew", pady=3)
        r += 1
        ttk.Label(settings, text="Y-axis label").grid(row=r, column=0, sticky="w", pady=3)
        ttk.Entry(settings, textvariable=self.ylabel, width=24).grid(row=r, column=1, sticky="ew", pady=3)
        r += 1
        ttk.Label(settings, text="Legend location").grid(row=r, column=0, sticky="w", pady=3)
        ttk.Combobox(settings, textvariable=self.legend_location, width=18, state="readonly", values=(
            "best", "upper right", "upper left", "lower left", "lower right", "right", "center left", "center right", "lower center", "upper center", "center"
        )).grid(row=r, column=1, sticky="w", pady=3)
        r += 1

        ttk.Separator(settings).grid(row=r, column=0, columnspan=2, sticky="ew", pady=10)
        r += 1
        ttk.Button(settings, text="Plot interactive", command=self.plot_interactive).grid(row=r, column=0, columnspan=2, sticky="ew", pady=(0, 6))
        r += 1
        ttk.Button(settings, text="Save figure…", command=self.save_figure).grid(row=r, column=0, columnspan=2, sticky="ew", pady=(0, 6))
        r += 1
        ttk.Label(settings, text="Interactive plot: click legend entries to hide/show curves; left-click an energy to print interpolated values; press 'a' to show all; 'n' to clear clicked points.", wraplength=260, justify="left").grid(row=r, column=0, columnspan=2, sticky="ew", pady=(8, 0))

        status = ttk.Frame(root)
        status.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(8, 0))
        status.columnconfigure(0, weight=1)
        ttk.Label(status, textvariable=self.status).grid(row=0, column=0, sticky="w")

    def _selected_indices(self) -> list[int]:
        return sorted([int(iid) for iid in self.tree.selection()])

    def _on_select(self, event=None):
        indices = self._selected_indices()
        if not indices:
            self.selected_file.set("No spectrum selected")
            self.selected_label.set("")
            return
        idx = indices[0]
        self.selected_file.set(self.files[idx].name)
        self.selected_label.set(self.labels[idx])

    def _focus_label_entry(self):
        self.label_entry.focus_set()
        self.label_entry.selection_range(0, "end")

    def add_files(self):
        paths = filedialog.askopenfilenames(
            title="Select processed ASTRA .dat files",
            filetypes=[("Processed spectra", "*.dat *.txt"), ("All files", "*.*")],
        )
        for p in paths:
            path = Path(p)
            if path in self.files:
                continue
            self.files.append(path)
            self.labels.append(self._default_label(path))
        self._refresh_tree(select_last=True)

    @staticmethod
    def _default_label(path: Path) -> str:
        return path.stem.replace("_norm", "").replace("_raw", "").replace("_flat", "")

    def remove_selected(self):
        selected = set(self._selected_indices())
        self.files = [f for i, f in enumerate(self.files) if i not in selected]
        self.labels = [lab for i, lab in enumerate(self.labels) if i not in selected]
        self._refresh_tree()

    def clear_files(self):
        self.files.clear()
        self.labels.clear()
        self._refresh_tree()

    def auto_labels(self):
        self.labels = [self._default_label(p) for p in self.files]
        self._refresh_tree(keep_selection=True)

    def apply_label(self):
        indices = self._selected_indices()
        if not indices:
            messagebox.showinfo("No selection", "Select one spectrum first.")
            return
        idx = indices[0]
        text = self.selected_label.get().strip()
        self.labels[idx] = text or self._default_label(self.files[idx])
        self._refresh_tree(select_indices=indices)

    def move_up(self):
        indices = self._selected_indices()
        if len(indices) != 1:
            messagebox.showinfo("Move up", "Select exactly one row to move.")
            return
        i = indices[0]
        if i == 0:
            return
        self.files[i - 1], self.files[i] = self.files[i], self.files[i - 1]
        self.labels[i - 1], self.labels[i] = self.labels[i], self.labels[i - 1]
        self._refresh_tree(select_indices=[i - 1])

    def move_down(self):
        indices = self._selected_indices()
        if len(indices) != 1:
            messagebox.showinfo("Move down", "Select exactly one row to move.")
            return
        i = indices[0]
        if i >= len(self.files) - 1:
            return
        self.files[i + 1], self.files[i] = self.files[i], self.files[i + 1]
        self.labels[i + 1], self.labels[i] = self.labels[i], self.labels[i + 1]
        self._refresh_tree(select_indices=[i + 1])

    def sort_by_filename(self):
        combined = sorted(zip(self.files, self.labels), key=lambda x: x[0].name.lower())
        self.files = [x[0] for x in combined]
        self.labels = [x[1] for x in combined]
        self._refresh_tree()
        
    def select_all(self):
        self.tree.selection_set(self.tree.get_children())
        self._on_select()
        self.status.set(f"{len(self.files)} file(s) selected for plotting.")

    def deselect_all(self):
        self.tree.selection_remove(self.tree.selection())
        self._on_select()
        self.status.set("No spectra selected for plotting.")

    def invert_selection(self):
        all_items = set(self.tree.get_children())
        selected = set(self.tree.selection())
        inverted = sorted(all_items - selected, key=int)
        self.tree.selection_set(inverted)
        self._on_select()
        self.status.set(f"{len(inverted)} file(s) selected for plotting.")

    def _refresh_tree(self, select_last: bool = False, keep_selection: bool = False, select_indices: list[int] | None = None):
        old_selection = self._selected_indices() if keep_selection else []
        self.tree.delete(*self.tree.get_children())
        for i, (p, lab) in enumerate(zip(self.files, self.labels)):
            self.tree.insert("", "end", iid=str(i), values=(p.name, lab))
        if select_indices is None:
            select_indices = old_selection
        if select_last and self.files:
            select_indices = [len(self.files) - 1]
        valid = [str(i) for i in select_indices if 0 <= i < len(self.files)]
        if valid:
            self.tree.selection_set(valid)
            self.tree.see(valid[0])
        self._on_select()
        self.status.set(f"{len(self.files)} file(s) loaded.")

    def save_session(self):
        path = filedialog.asksaveasfilename(
            title="Save viewer session",
            defaultextension=".json",
            filetypes=[("JSON", "*.json"), ("All files", "*.*")],
        )
        if not path:
            return
        payload = {
            "files": [str(p) for p in self.files],
            "labels": self.labels,
            "settings": {
                "e_min": self.e_min.get(), "e_max": self.e_max.get(),
                "plot_channel": self.plot_channel.get(),
                "use_smoothing": self.use_smoothing.get(), "show_raw_too": self.show_raw_too.get(),
                "sg_window": self.sg_window.get(), "sg_polyorder": self.sg_polyorder.get(),
                "linewidth": self.linewidth.get(), "legend_fontsize": self.legend_fontsize.get(),
                "fig_width": self.fig_width.get(), "fig_height": self.fig_height.get(),
                "xlabel": self.xlabel.get(), "ylabel": self.ylabel.get(),
                "legend_location": self.legend_location.get(), "title_text": self.title_text.get(),
                "grid_on": self.grid_on.get(),
            },
        }
        Path(path).write_text(json.dumps(payload, indent=2), encoding="utf-8")
        self.status.set(f"Saved session: {path}")

    def load_session(self):
        path = filedialog.askopenfilename(
            title="Load viewer session",
            filetypes=[("JSON", "*.json"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            payload = json.loads(Path(path).read_text(encoding="utf-8"))
            self.files = [Path(p) for p in payload.get("files", [])]
            self.labels = list(payload.get("labels", [self._default_label(p) for p in self.files]))
            settings = payload.get("settings", {})
            for key, var in [
                ("e_min", self.e_min), ("e_max", self.e_max), ("plot_channel", self.plot_channel), ("sg_window", self.sg_window),
                ("sg_polyorder", self.sg_polyorder), ("linewidth", self.linewidth),
                ("legend_fontsize", self.legend_fontsize), ("fig_width", self.fig_width),
                ("fig_height", self.fig_height), ("xlabel", self.xlabel), ("ylabel", self.ylabel),
                ("legend_location", self.legend_location), ("title_text", self.title_text),
            ]:
                if key in settings:
                    var.set(settings[key])
            self.use_smoothing.set(bool(settings.get("use_smoothing", self.use_smoothing.get())))
            self.show_raw_too.set(bool(settings.get("show_raw_too", self.show_raw_too.get())))
            self.grid_on.set(bool(settings.get("grid_on", self.grid_on.get())))
            self._refresh_tree()
        except Exception as exc:
            messagebox.showerror("Could not load session", str(exc))

    def _float(self, name: str, var: tk.StringVar) -> float:
        try:
            return float(var.get())
        except ValueError as exc:
            raise ValueError(f"Invalid {name}: {var.get()!r}") from exc

    def _int(self, name: str, var: tk.StringVar) -> int:
        try:
            return int(var.get())
        except ValueError as exc:
            raise ValueError(f"Invalid {name}: {var.get()!r}") from exc

    def _load_spectra(self):
        if not self.files:
            raise ValueError("No files selected.")

        selected = self._selected_indices()
        if not selected:
            selected = list(range(len(self.files)))

        e_min = self._float("energy minimum", self.e_min)
        e_max = self._float("energy maximum", self.e_max)

        if e_min >= e_max:
            raise ValueError("Energy minimum must be smaller than energy maximum.")

        sg_window = self._int("SG window", self.sg_window)
        sg_polyorder = self._int("SG polyorder", self.sg_polyorder)

        spectra = []
        for i in selected:
            path = self.files[i]
            label = self.labels[i]

            energy, y, y_col = read_spectrum_dat(path, channel=self.plot_channel.get())

            mask = (energy >= e_min) & (energy <= e_max)
            x_plot = energy[mask]
            y_raw = y[mask]

            if len(x_plot) == 0:
                raise ValueError(f"No points in selected energy range for {path.name}")

            if self.use_smoothing.get():
                y_plot = _safe_savgol(y_raw, sg_window, sg_polyorder)
            else:
                y_plot = y_raw.copy()

            spectra.append({
                "path": path,
                "label": label,
                "y_col": y_col,
                "energy": x_plot,
                "y_raw": y_raw,
                "y_plot": y_plot
            })

        return spectra
    
    def _make_figure(self, interactive: bool = True):
        import matplotlib
        try:
            matplotlib.use("TkAgg")
        except Exception:
            pass
        import matplotlib.pyplot as plt

        spectra = self._load_spectra()
        fig_w = self._float("figure width", self.fig_width)
        fig_h = self._float("figure height", self.fig_height)
        lw = self._float("line width", self.linewidth)
        legend_fs = self._float("legend font size", self.legend_fontsize)

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        lines = []
        raw_lines = []
        for spec in spectra:
            if self.show_raw_too.get() and self.use_smoothing.get():
                raw_line, = ax.plot(spec["energy"], spec["y_raw"], lw=1.2, ls="--", alpha=0.45, label=f"{spec['label']} (raw)")
                raw_lines.append(raw_line)
            line, = ax.plot(spec["energy"], spec["y_plot"], lw=lw, label=spec["label"])
            lines.append(line)

        ax.set_xlim(self._float("energy minimum", self.e_min), self._float("energy maximum", self.e_max))
        ax.set_xlabel(self.xlabel.get().strip() or "Energy (eV)", fontsize=13)
        ax.set_ylabel(self.ylabel.get().strip() or "Signal", fontsize=13)
        title = self.title_text.get().strip()
        if title:
            ax.set_title(title, fontsize=15)
        if self.grid_on.get():
            ax.grid(True, alpha=0.25)
        for spine in ax.spines.values():
            spine.set_linewidth(2.0)
        ax.tick_params(axis="both", which="major", direction="out", length=6, width=1.8, labelsize=12)
        ax.minorticks_on()
        ax.tick_params(axis="both", which="minor", direction="out", length=3, width=1.2)
        legend = ax.legend(loc=self.legend_location.get(), frameon=False, fontsize=legend_fs)

        all_lines = lines + raw_lines
        if interactive:
            self._attach_interactions(fig, ax, legend, all_lines)
        fig.tight_layout()
        return fig

    def _attach_interactions(self, fig, ax, legend, all_lines):
        import numpy as np
        import matplotlib.pyplot as plt

        lined = {}
        for legline, origline in zip(legend.get_lines(), all_lines):
            legline.set_picker(5)
            lined[legline] = origline

        annot = ax.annotate("", xy=(0, 0), xytext=(15, 15), textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="white", ec="black", alpha=0.85),
                            arrowprops=dict(arrowstyle="->"))
        annot.set_visible(False)
        clicked_points = []
        scatter = ax.scatter([], [], c="red", s=35, zorder=6)

        def on_pick(event):
            legline = event.artist
            if legline not in lined:
                return
            origline = lined[legline]
            visible = not origline.get_visible()
            origline.set_visible(visible)
            legline.set_alpha(1.0 if visible else 0.2)
            for i, lgl in enumerate(legend.get_lines()):
                if lgl is legline:
                    legend.get_texts()[i].set_alpha(1.0 if visible else 0.2)
                    break
            fig.canvas.draw_idle()

        def update_annot(line, idx):
            xdata, ydata = line.get_data()
            annot.xy = (xdata[idx], ydata[idx])
            annot.set_text(f"{line.get_label()}\nE = {xdata[idx]:.3f} eV\nY = {ydata[idx]:.6f}")

        def on_move(event):
            if event.inaxes != ax or event.xdata is None or event.ydata is None:
                if annot.get_visible():
                    annot.set_visible(False)
                    fig.canvas.draw_idle()
                return
            nearest_line = None
            nearest_idx = None
            min_dist = np.inf
            y0, y1 = ax.get_ylim(); x0, x1 = ax.get_xlim()
            y_scale = max(abs(y1 - y0), 1e-12); x_scale = max(abs(x1 - x0), 1e-12)
            for line in all_lines:
                if not line.get_visible():
                    continue
                xdata = line.get_xdata(); ydata = line.get_ydata()
                if len(xdata) == 0:
                    continue
                idx = int(np.argmin(np.abs(xdata - event.xdata)))
                dx = (xdata[idx] - event.xdata) / x_scale
                dy = (ydata[idx] - event.ydata) / y_scale
                dist = dx * dx + dy * dy
                if dist < min_dist:
                    min_dist = dist; nearest_line = line; nearest_idx = idx
            if nearest_line is not None and min_dist < 0.002:
                update_annot(nearest_line, nearest_idx)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            elif annot.get_visible():
                annot.set_visible(False)
                fig.canvas.draw_idle()

        def on_click(event):
            if event.inaxes != ax or event.button != 1 or event.xdata is None:
                return
            tb = getattr(fig.canvas, "toolbar", None)
            if tb is not None and getattr(tb, "mode", "") != "":
                return
            e_click = float(event.xdata)
            print(f"\nE = {e_click:.3f} eV")
            new_points = []
            for line in all_lines:
                if not line.get_visible():
                    continue
                xdata = np.asarray(line.get_xdata(), dtype=float)
                ydata = np.asarray(line.get_ydata(), dtype=float)
                if len(xdata) == 0:
                    continue
                y_val = float(np.interp(e_click, xdata, ydata))
                print(f"{line.get_label()}: {y_val:.6f}")
                new_points.append([e_click, y_val])
            if new_points:
                clicked_points.extend(new_points)
                scatter.set_offsets(np.asarray(clicked_points))
                fig.canvas.draw_idle()

        def on_key(event):
            if event.key in ("enter", "escape"):
                plt.close(fig)
            elif event.key == "a":
                for line in all_lines:
                    line.set_visible(True)
                for legline in legend.get_lines():
                    legline.set_alpha(1.0)
                for text in legend.get_texts():
                    text.set_alpha(1.0)
                fig.canvas.draw_idle()
            elif event.key == "n":
                clicked_points.clear()
                scatter.set_offsets(np.empty((0, 2)))
                fig.canvas.draw_idle()

        fig.canvas.mpl_connect("pick_event", on_pick)
        fig.canvas.mpl_connect("motion_notify_event", on_move)
        fig.canvas.mpl_connect("button_press_event", on_click)
        fig.canvas.mpl_connect("key_press_event", on_key)

    def plot_interactive(self):
        try:
            import matplotlib.pyplot as plt
            fig = self._make_figure(interactive=True)
            plt.show(block=False)
            self.status.set("Interactive plot opened.")
        except Exception as exc:
            messagebox.showerror("Could not plot", str(exc))

    def save_figure(self):
        if not self.files:
            messagebox.showerror("No files", "Select spectra first.")
            return
        path = filedialog.asksaveasfilename(
            title="Save spectrum comparison plot",
            defaultextension=".png",
            filetypes=[("PNG image", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            import matplotlib.pyplot as plt
            fig = self._make_figure(interactive=False)
            fig.savefig(path, dpi=300)
            plt.close(fig)
            self.status.set(f"Saved: {path}")
        except Exception as exc:
            messagebox.showerror("Could not save plot", str(exc))


def open_spectrum_viewer(master=None):
    viewer = SpectrumViewer(master)
    viewer.focus_set()
    return viewer


def main():
    app = tk.Tk()
    app.withdraw()
    viewer = SpectrumViewer(app)
    viewer.protocol("WM_DELETE_WINDOW", app.destroy)
    app.mainloop()


if __name__ == "__main__":
    main()
