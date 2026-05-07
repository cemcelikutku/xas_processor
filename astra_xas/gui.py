from __future__ import annotations

import json
import time
import queue
import threading
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox
from tkinter import ttk

from .config import AstraConfig
from .processor import process_folder
from .spectrum_viewer import open_spectrum_viewer


class AstraGui(tk.Tk):
    """Tkinter interface for ASTRA XAS Processor."""

    def __init__(self):
        super().__init__()

        self.title("ASTRA XAS Processor")
        self.geometry("1120x780")
        self.minsize(980, 640)

        self.input_dir = tk.StringVar()
        self.output_dir = tk.StringVar()
        self.analysis_mode = tk.StringVar(value="fluo")
        self.alignment_source = tk.StringVar(value="inline_ref")
        self.foil_mode = tk.StringVar(value="trans")
        self.foil_keyword = tk.StringVar(value="foil")
        self.e0 = tk.StringVar(value="7121.030")
        self.pre1 = tk.StringVar(value="-229.740")
        self.pre2 = tk.StringVar(value="-49.980")
        self.norm1 = tk.StringVar(value="55.070")
        self.norm2 = tk.StringVar(value="227.220")
        self.nnorm = tk.StringVar(value="1")
        self.align_min = tk.StringVar(value="7100.0")
        self.align_max = tk.StringVar(value="7140.0")
        self.shift_min = tk.StringVar(value="-5.0")
        self.shift_max = tk.StringVar(value="5.0")
        self.alignment_quality_warn_threshold = tk.StringVar(value="0.7")
        self.alignment_grid_points = tk.StringVar(value="50")
        self.save_drift_plot = tk.BooleanVar(value=False)
        self.fluo_factor = tk.StringVar(value="1e-11")
        self.enable_deglitching = tk.BooleanVar(value=False)
        self.deglitch_mode = tk.StringVar(value="automatic")
        self.enable_auto_deglitch = tk.BooleanVar(value=False)
        self.deglitch_method = tk.StringVar(value="interpolate")
        self.deglitch_threshold = tk.StringVar(value="5.0")
        self.deglitch_window = tk.StringVar(value="5")
        self.deglitch_min_energy = tk.StringVar(value="")
        self.deglitch_max_energy = tk.StringVar(value="")

        self.enable_manual_deglitch_range = tk.BooleanVar(value=False)
        self.manual_deglitch_min_energy = tk.StringVar(value="")
        self.manual_deglitch_max_energy = tk.StringVar(value="")
        self.manual_deglitch_margin_points = tk.StringVar(value="5")

        self.plot_detector_raw_overview = tk.BooleanVar(value=False)
        self.plot_processed_overview = tk.BooleanVar(value=True)
        self.plot_bkgcorr_overview = tk.BooleanVar(value=False)
        self.plot_norm_overview = tk.BooleanVar(value=True)
        self.plot_replicate_qc = tk.BooleanVar(value=True)
        # Backward-compatible internal alias. Do not expose as "raw" in the GUI.
        self.plot_raw_overview = self.plot_processed_overview
        self.plot_min = tk.StringVar(value="7100.0")
        self.plot_max = tk.StringVar(value="7160.0")

        self.status = tk.StringVar(value="Ready")
        self._log_queue: queue.Queue[str] = queue.Queue()
        self._running = False
        self._suppress_log_until = 0
        self._auto_deglitch_widgets = []
        self._manual_deglitch_widgets = []
        self._deglitch_mode_widgets = []
        self._build()
        self.alignment_source.trace_add("write", self._update_alignment_ui)
        self.enable_deglitching.trace_add("write", self._update_deglitch_ui)
        self.deglitch_mode.trace_add("write", self._update_deglitch_ui)
        self._update_alignment_ui()
        self._update_deglitch_ui()
        self.after(100, self._drain_log_queue)

    def _build(self):
        style = ttk.Style(self)
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass

        style.configure("Title.TLabel", font=("TkDefaultFont", 16, "bold"))
        style.configure("Subtitle.TLabel", font=("TkDefaultFont", 10))
        style.configure("Section.TLabelframe.Label", font=("TkDefaultFont", 10, "bold"))
        style.configure("Run.TButton", font=("TkDefaultFont", 10, "bold"), padding=8)

        container = ttk.Frame(self)
        container.pack(fill="both", expand=True)

        canvas = tk.Canvas(container, highlightthickness=0)
        scrollbar = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar.set)

        scrollbar.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        scrollable = ttk.Frame(canvas)
        window_id = canvas.create_window((0, 0), window=scrollable, anchor="nw")

        def on_frame_configure(_event=None):
            canvas.configure(scrollregion=canvas.bbox("all"))

        def on_canvas_configure(event):
            canvas.itemconfigure(window_id, width=event.width)

        scrollable.bind("<Configure>", on_frame_configure)
        canvas.bind("<Configure>", on_canvas_configure)

        def on_mousewheel(event):
            canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        canvas.bind_all("<MouseWheel>", on_mousewheel)
        canvas.bind_all("<Button-4>", lambda _e: canvas.yview_scroll(-1, "units"))
        canvas.bind_all("<Button-5>", lambda _e: canvas.yview_scroll(1, "units"))

        root = ttk.Frame(scrollable, padding=14)
        root.pack(fill="both", expand=True)

        root.columnconfigure(0, weight=0)
        root.columnconfigure(1, weight=1)
        root.rowconfigure(1, weight=1)

        header = ttk.Frame(root)
        header.grid(row=0, column=0, columnspan=2, sticky="ew", pady=(0, 12))
        header.columnconfigure(0, weight=1)

        ttk.Label(
            header,
            text="🦇 ASTRA XAS Processor",
            style="Title.TLabel",
        ).grid(row=0, column=0, sticky="w")

        ttk.Label(
            header,
            text="Foil drift correction, replicate averaging, Athena-like normalization, and spectrum visualization for ASTRA .xasd data.",
            style="Subtitle.TLabel",
        ).grid(row=1, column=0, sticky="w", pady=(3, 0))

        left = ttk.Frame(root)
        left.grid(row=1, column=0, sticky="nsew", padx=(0, 14))

        right = ttk.Frame(root)
        right.grid(row=1, column=1, sticky="nsew")
        right.rowconfigure(1, weight=1)
        right.columnconfigure(0, weight=1)

        self._build_input_section(left)
        self._build_basic_section(left)
        self._build_advanced_section(left)
        self._build_deglitch_section(left)
        self._build_plot_section(left)
        self._build_buttons(left)
        self._build_log_section(right)

        statusbar = ttk.Frame(root)
        statusbar.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(10, 0))
        statusbar.columnconfigure(0, weight=1)

        ttk.Label(statusbar, textvariable=self.status).grid(row=0, column=0, sticky="w")
        self.progress = ttk.Progressbar(statusbar, mode="indeterminate", length=180)
        self.progress.grid(row=0, column=1, sticky="e")

    def _build_input_section(self, parent):
        frame = ttk.LabelFrame(parent, text="1. Folders", padding=10, style="Section.TLabelframe")
        frame.pack(fill="x", pady=(0, 10))
        frame.columnconfigure(1, weight=1)

        ttk.Label(frame, text="Input folder").grid(row=0, column=0, sticky="w", pady=4)
        ttk.Entry(frame, textvariable=self.input_dir, width=52).grid(
            row=0, column=1, sticky="ew", padx=6, pady=4
        )
        ttk.Button(frame, text="Browse…", command=self.pick_input).grid(row=0, column=2, pady=4)

        ttk.Label(frame, text="Output folder").grid(row=1, column=0, sticky="w", pady=4)
        ttk.Entry(frame, textvariable=self.output_dir, width=52).grid(
            row=1, column=1, sticky="ew", padx=6, pady=4
        )
        ttk.Button(frame, text="Browse…", command=self.pick_output).grid(row=1, column=2, pady=4)

        ttk.Label(
            frame,
            text="Leave output empty to create <input-folder>-processed automatically.",
        ).grid(row=2, column=1, columnspan=2, sticky="w", padx=6, pady=(0, 2))

    def _build_basic_section(self, parent):
        frame = ttk.LabelFrame(parent, text="2. Main settings", padding=10, style="Section.TLabelframe")
        frame.pack(fill="x", pady=(0, 10))
        frame.columnconfigure(1, weight=1)

        ttk.Label(frame, text="Sample signal").grid(row=0, column=0, sticky="w", pady=4)
        ttk.Combobox(
            frame,
            textvariable=self.analysis_mode,
            values=("fluo", "trans", "ref"),
            width=12,
            state="readonly",
        ).grid(row=0, column=1, sticky="w", padx=6, pady=4)

        ttk.Label(frame, text="Alignment source").grid(row=1, column=0, sticky="w", pady=4)
        ttk.Combobox(
            frame,
            textvariable=self.alignment_source,
            values=("inline_ref", "separate_foil"),
            width=16,
            state="readonly",
        ).grid(row=1, column=1, sticky="w", padx=6, pady=4)

        ttk.Label(frame, text="Foil alignment signal").grid(row=2, column=0, sticky="w", pady=4)
        ttk.Combobox(
            frame,
            textvariable=self.foil_mode,
            values=("trans", "ref", "fluo"),
            width=12,
            state="readonly",
        ).grid(row=2, column=1, sticky="w", padx=6, pady=4)

        ttk.Label(frame, text="Foil filename keyword").grid(row=3, column=0, sticky="w", pady=4)
        self.foil_keyword_entry = ttk.Entry(frame, textvariable=self.foil_keyword, width=16)
        self.foil_keyword_entry.grid(row=3, column=1, sticky="w", padx=6, pady=4)

        ttk.Label(frame, text="E0 / eV").grid(row=4, column=0, sticky="w", pady=4)
        ttk.Entry(frame, textvariable=self.e0, width=16).grid(
            row=4, column=1, sticky="w", padx=6, pady=4
        )

        ttk.Label(
            frame,
            text="Parameters are user-defined. Save a config file for each edge or experiment type.",
            wraplength=430,
        ).grid(row=5, column=0, columnspan=2, sticky="w", pady=(6, 0))

    def _build_advanced_section(self, parent):
        frame = ttk.LabelFrame(
            parent,
            text="3. Advanced normalization / alignment",
            padding=10,
            style="Section.TLabelframe",
        )
        frame.pack(fill="x", pady=(0, 10))
        frame.columnconfigure(1, weight=1)
        frame.columnconfigure(3, weight=1)

        rows = [
            ("pre1", self.pre1, "pre2", self.pre2),
            ("norm1", self.norm1, "norm2", self.norm2),
            ("normalization order", self.nnorm, "", None),
            ("align min / eV", self.align_min, "align max / eV", self.align_max),
            ("shift min / eV", self.shift_min, "shift max / eV", self.shift_max),
            ("Quality warn threshold (0–1)", self.alignment_quality_warn_threshold, "Grid points (coarse search)", self.alignment_grid_points),
            ("fluo factor", self.fluo_factor, "", None),
        ]

        for r, (label1, var1, label2, var2) in enumerate(rows):
            ttk.Label(frame, text=label1).grid(row=r, column=0, sticky="w", pady=3)
            ttk.Entry(frame, textvariable=var1, width=14).grid(
                row=r, column=1, sticky="w", padx=6, pady=3
            )

            if var2 is not None:
                ttk.Label(frame, text=label2).grid(row=r, column=2, sticky="w", padx=(10, 0), pady=3)
                ttk.Entry(frame, textvariable=var2, width=14).grid(
                    row=r, column=3, sticky="w", padx=6, pady=3
                )

        ttk.Checkbutton(
            frame,
            text="Save drift plot",
            variable=self.save_drift_plot,
        ).grid(row=len(rows), column=0, columnspan=4, sticky="w", pady=(5, 0))

    def _build_deglitch_section(self, parent):
        frame = ttk.LabelFrame(
            parent,
            text="3a. Deglitching",
            padding=10,
            style="Section.TLabelframe",
        )
        frame.pack(fill="x", pady=(0, 10))
        frame.columnconfigure(1, weight=1)
        frame.columnconfigure(3, weight=1)

        ttk.Checkbutton(
            frame,
            text="Enable deglitching",
            variable=self.enable_deglitching,
        ).grid(row=0, column=0, columnspan=4, sticky="w", pady=(0, 8))

        ttk.Label(frame, text="Mode").grid(row=1, column=0, sticky="w", pady=3)
        mode_combo = ttk.Combobox(
            frame,
            textvariable=self.deglitch_mode,
            values=("automatic", "manual", "both"),
            width=14,
            state="readonly",
        )
        mode_combo.grid(row=1, column=1, sticky="w", padx=6, pady=3)
        self._deglitch_mode_widgets.append(mode_combo)

        auto_label = ttk.Label(frame, text="Automatic point spikes")
        auto_label.grid(row=2, column=0, columnspan=4, sticky="w", pady=(8, 3))
        self._auto_deglitch_widgets.append(auto_label)

        threshold_label = ttk.Label(frame, text="Threshold")
        threshold_label.grid(row=3, column=0, sticky="w", pady=3)
        threshold_entry = ttk.Entry(frame, textvariable=self.deglitch_threshold, width=14)
        threshold_entry.grid(row=3, column=1, sticky="w", padx=6, pady=3)
        self._auto_deglitch_widgets.extend([threshold_label, threshold_entry])

        window_label = ttk.Label(frame, text="Window half-width")
        window_label.grid(row=4, column=0, sticky="w", pady=3)
        window_entry = ttk.Entry(frame, textvariable=self.deglitch_window, width=14)
        window_entry.grid(
            row=4, column=1, sticky="w", padx=6, pady=3
        )
        self._auto_deglitch_widgets.extend([window_label, window_entry])

        auto_min_label = ttk.Label(frame, text="Min energy / eV")
        auto_min_label.grid(row=4, column=2, sticky="w", padx=(10, 0), pady=3)
        auto_min_entry = ttk.Entry(frame, textvariable=self.deglitch_min_energy, width=14)
        auto_min_entry.grid(
            row=4, column=3, sticky="w", padx=6, pady=3
        )
        self._auto_deglitch_widgets.extend([auto_min_label, auto_min_entry])

        auto_max_label = ttk.Label(frame, text="Max energy / eV")
        auto_max_label.grid(row=5, column=0, sticky="w", pady=3)
        auto_max_entry = ttk.Entry(frame, textvariable=self.deglitch_max_energy, width=14)
        auto_max_entry.grid(
            row=5, column=1, sticky="w", padx=6, pady=3
        )
        self._auto_deglitch_widgets.extend([auto_max_label, auto_max_entry])

        manual_label = ttk.Label(frame, text="Manual range")
        manual_label.grid(row=6, column=0, columnspan=4, sticky="w", pady=(8, 3))
        self._manual_deglitch_widgets.append(manual_label)

        manual_min_label = ttk.Label(frame, text="Min energy / eV")
        manual_min_label.grid(row=7, column=0, sticky="w", pady=3)
        manual_min_entry = ttk.Entry(frame, textvariable=self.manual_deglitch_min_energy, width=14)
        manual_min_entry.grid(
            row=7, column=1, sticky="w", padx=6, pady=3
        )
        self._manual_deglitch_widgets.extend([manual_min_label, manual_min_entry])

        manual_max_label = ttk.Label(frame, text="Max energy / eV")
        manual_max_label.grid(row=7, column=2, sticky="w", padx=(10, 0), pady=3)
        manual_max_entry = ttk.Entry(frame, textvariable=self.manual_deglitch_max_energy, width=14)
        manual_max_entry.grid(
            row=7, column=3, sticky="w", padx=6, pady=3
        )
        self._manual_deglitch_widgets.extend([manual_max_label, manual_max_entry])

        margin_label = ttk.Label(frame, text="Margin points")
        margin_label.grid(row=8, column=0, sticky="w", pady=3)
        margin_entry = ttk.Entry(frame, textvariable=self.manual_deglitch_margin_points, width=14)
        margin_entry.grid(
            row=8, column=1, sticky="w", padx=6, pady=3
        )
        self._manual_deglitch_widgets.extend([margin_label, margin_entry])

    def _build_plot_section(self, parent):
        frame = ttk.LabelFrame(parent, text="4. Automatic plots", padding=10, style="Section.TLabelframe")
        frame.pack(fill="x", pady=(0, 10))
        frame.columnconfigure(1, weight=1)

        ttk.Checkbutton(
            frame,
            text="Aligned averaged IF overview",
            variable=self.plot_detector_raw_overview,
        ).grid(row=0, column=0, columnspan=2, sticky="w", pady=2)

        ttk.Checkbutton(
            frame,
            text="Processed μ(E) overview",
            variable=self.plot_processed_overview,
        ).grid(row=1, column=0, columnspan=2, sticky="w", pady=2)

        ttk.Checkbutton(
            frame,
            text="Background-corrected overview",
            variable=self.plot_bkgcorr_overview,
        ).grid(row=2, column=0, columnspan=2, sticky="w", pady=2)

        ttk.Checkbutton(
            frame,
            text="Normalized overview",
            variable=self.plot_norm_overview,
        ).grid(row=3, column=0, columnspan=2, sticky="w", pady=2)

        ttk.Checkbutton(
            frame,
            text="Replicate QC plots",
            variable=self.plot_replicate_qc,
        ).grid(row=4, column=0, columnspan=2, sticky="w", pady=2)

        range_frame = ttk.Frame(frame)
        range_frame.grid(row=5, column=0, columnspan=2, sticky="w", pady=(6, 0))

        ttk.Label(range_frame, text="Plot energy range / eV").grid(row=0, column=0, sticky="w")
        ttk.Entry(range_frame, textvariable=self.plot_min, width=10).grid(row=0, column=1, padx=(8, 4))
        ttk.Label(range_frame, text="to").grid(row=0, column=2)
        ttk.Entry(range_frame, textvariable=self.plot_max, width=10).grid(row=0, column=3, padx=(4, 0))

    def _build_buttons(self, parent):
        frame = ttk.Frame(parent)
        frame.pack(fill="x", pady=(0, 10))
        frame.columnconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)
        frame.columnconfigure(2, weight=1)

        self.run_button = ttk.Button(
            frame,
            text="Run processing",
            command=self.run_processing,
            style="Run.TButton",
        )
        self.run_button.grid(row=0, column=0, columnspan=3, sticky="ew", pady=(0, 8))

        ttk.Button(frame, text="Load config", command=self.load_config).grid(
            row=1, column=0, sticky="ew", padx=(0, 4)
        )
        ttk.Button(frame, text="Save config", command=self.save_config).grid(
            row=1, column=1, sticky="ew", padx=4
        )
        ttk.Button(frame, text="Clear log", command=self.clear_log).grid(
            row=1, column=2, sticky="ew", padx=(4, 0)
        )

        ttk.Button(
            frame,
            text="Open Spectrum Viewer",
            command=lambda: open_spectrum_viewer(self),
        ).grid(row=2, column=0, columnspan=3, sticky="ew", pady=(8, 0))

    def _build_log_section(self, parent):
        top = ttk.Frame(parent)
        top.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        top.columnconfigure(0, weight=1)

        ttk.Label(top, text="Processing log", style="Section.TLabelframe.Label").grid(
            row=0, column=0, sticky="w"
        )

        log_frame = ttk.Frame(parent)
        log_frame.grid(row=1, column=0, sticky="nsew")
        log_frame.rowconfigure(0, weight=1)
        log_frame.columnconfigure(0, weight=1)

        self.log_box = tk.Text(
            log_frame,
            wrap="word",
            height=24,
            width=78,
            state="disabled",
        )
        self.log_box.grid(row=0, column=0, sticky="nsew")

        scroll = ttk.Scrollbar(log_frame, orient="vertical", command=self.log_box.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        self.log_box.configure(yscrollcommand=scroll.set)

        help_box = ttk.LabelFrame(parent, text="Expected workflow", padding=10)
        help_box.grid(row=2, column=0, sticky="ew", pady=(10, 0))

        ttk.Label(
            help_box,
            text=(
                "1) Select the folder containing ASTRA .xasd files.  "
                "2) Keep sample signal as fluo for fluorescence XAS.  "
                "3) Choose inline_ref if I1/I2 is measured in each scan, or separate_foil for separate foil files.  "
                "4) Run and check ASTRA_processing_report.txt and the plots folder."
            ),
            wraplength=520,
            justify="left",
        ).pack(fill="x")

    def _update_alignment_ui(self, *args):
        """Disable separate-foil-only controls when inline reference alignment is selected."""
        if hasattr(self, "foil_keyword_entry"):
            if self.alignment_source.get() == "inline_ref":
                self.foil_keyword_entry.configure(state="disabled")
            else:
                self.foil_keyword_entry.configure(state="normal")

    def _deglitch_enable_flags(self) -> tuple[bool, bool]:
        if not self.enable_deglitching.get():
            return False, False
        mode = self.deglitch_mode.get()
        if mode not in {"automatic", "manual", "both"}:
            mode = "automatic"
        return mode in {"automatic", "both"}, mode in {"manual", "both"}

    def _sync_deglitch_backend_vars(self):
        auto_enabled, manual_enabled = self._deglitch_enable_flags()
        self.enable_auto_deglitch.set(auto_enabled)
        self.enable_manual_deglitch_range.set(manual_enabled)

    def _set_widget_state(self, widgets, state: str):
        for widget in widgets:
            try:
                widget.configure(state=state)
            except tk.TclError:
                pass

    def _update_deglitch_ui(self, *args):
        master_enabled = self.enable_deglitching.get()
        auto_enabled, manual_enabled = self._deglitch_enable_flags()

        self._set_widget_state(
            self._deglitch_mode_widgets,
            "readonly" if master_enabled else "disabled",
        )
        self._set_widget_state(
            self._auto_deglitch_widgets,
            "normal" if auto_enabled else "disabled",
        )
        self._set_widget_state(
            self._manual_deglitch_widgets,
            "normal" if manual_enabled else "disabled",
        )
        self._sync_deglitch_backend_vars()

    def pick_input(self):
        d = filedialog.askdirectory(title="Select folder containing .xasd files")
        if d:
            self.input_dir.set(d)
            if not self.output_dir.get().strip():
                p = Path(d).expanduser().resolve()
                self.output_dir.set(str(p.parent / f"{p.name}-processed"))

    def pick_output(self):
        d = filedialog.askdirectory(title="Select output folder")
        if d:
            self.output_dir.set(d)

    def _float(self, name: str, var: tk.StringVar) -> float:
        try:
            return float(var.get())
        except ValueError as exc:
            raise ValueError(f"Invalid numeric value for {name}: {var.get()!r}") from exc

    def _float_or_none(self, name: str, var: tk.StringVar) -> float | None:
        value = var.get().strip()
        if value == "":
            return None
        try:
            return float(value)
        except ValueError as exc:
            raise ValueError(f"Invalid numeric value for {name}: {var.get()!r}") from exc

    def _int(self, name: str, var: tk.StringVar) -> int:
        try:
            return int(var.get())
        except ValueError as exc:
            raise ValueError(f"Invalid integer value for {name}: {var.get()!r}") from exc

    def build_config(self) -> AstraConfig:
        keyword = self.foil_keyword.get().strip()
        if self.alignment_source.get() == "separate_foil" and not keyword:
            raise ValueError("Foil filename keyword cannot be empty when alignment source is separate_foil.")

        align_min = self._float("align min", self.align_min)
        align_max = self._float("align max", self.align_max)
        shift_min = self._float("shift min", self.shift_min)
        shift_max = self._float("shift max", self.shift_max)
        plot_min = self._float("plot min", self.plot_min)
        plot_max = self._float("plot max", self.plot_max)
        nnorm = self._int("normalization order", self.nnorm)
        deglitch_threshold = self._float("deglitch threshold", self.deglitch_threshold)
        deglitch_window = self._int("deglitch window", self.deglitch_window)
        deglitch_min_energy = self._float_or_none("deglitch min energy", self.deglitch_min_energy)
        deglitch_max_energy = self._float_or_none("deglitch max energy", self.deglitch_max_energy)
        manual_deglitch_min_energy = self._float_or_none("manual deglitch min energy", self.manual_deglitch_min_energy)
        manual_deglitch_max_energy = self._float_or_none("manual deglitch max energy", self.manual_deglitch_max_energy)
        manual_deglitch_margin_points = self._int("manual deglitch margin points", self.manual_deglitch_margin_points)
        alignment_quality_warn_threshold = self._float(
            "alignment quality warn threshold",
            self.alignment_quality_warn_threshold,
        )
        alignment_grid_points = self._int(
            "alignment grid points",
            self.alignment_grid_points,
        )

        if align_min >= align_max:
            raise ValueError("align min must be smaller than align max.")
        if shift_min >= shift_max:
            raise ValueError("shift min must be smaller than shift max.")
        if plot_min >= plot_max:
            raise ValueError("plot min must be smaller than plot max.")
        if deglitch_min_energy is not None and deglitch_max_energy is not None and deglitch_min_energy >= deglitch_max_energy:
            raise ValueError("deglitch min energy must be smaller than deglitch max energy.")
        if manual_deglitch_min_energy is not None and manual_deglitch_max_energy is not None and manual_deglitch_min_energy >= manual_deglitch_max_energy:
            raise ValueError("manual deglitch min energy must be smaller than manual deglitch max energy.")
        if nnorm < 0:
            raise ValueError("normalization order must be 0 or a positive integer.")
        if deglitch_window < 1:
            raise ValueError("deglitch window must be a positive integer.")
        if manual_deglitch_margin_points < 1:
            raise ValueError("manual deglitch margin points must be a positive integer.")
        if alignment_grid_points < 5:
            raise ValueError("Alignment grid points must be at least 5.")
        if not (0.0 <= alignment_quality_warn_threshold <= 1.0):
            raise ValueError("Alignment quality warn threshold must be between 0.0 and 1.0.")

        auto_deglitch_enabled, manual_deglitch_enabled = self._deglitch_enable_flags()
        self.enable_auto_deglitch.set(auto_deglitch_enabled)
        self.enable_manual_deglitch_range.set(manual_deglitch_enabled)

        return AstraConfig(
            analysis_mode=self.analysis_mode.get(),
            alignment_source=self.alignment_source.get(),
            foil_alignment_mode=self.foil_mode.get(),
            foil_keyword=keyword,
            e0=self._float("E0", self.e0),
            pre1=self._float("pre1", self.pre1),
            pre2=self._float("pre2", self.pre2),
            norm1=self._float("norm1", self.norm1),
            norm2=self._float("norm2", self.norm2),
            nnorm=nnorm,
            align_window_min=align_min,
            align_window_max=align_max,
            shift_bound_min=shift_min,
            shift_bound_max=shift_max,
            alignment_quality_warn_threshold=alignment_quality_warn_threshold,
            alignment_grid_points=alignment_grid_points,
            fluo_multiplicative_constant=self._float("fluo factor", self.fluo_factor),
            enable_auto_deglitch=auto_deglitch_enabled,
            deglitch_method="interpolate",
            deglitch_threshold=deglitch_threshold,
            deglitch_window=deglitch_window,
            deglitch_min_energy=deglitch_min_energy,
            deglitch_max_energy=deglitch_max_energy,
            enable_manual_deglitch_range=manual_deglitch_enabled,
            manual_deglitch_min_energy=manual_deglitch_min_energy,
            manual_deglitch_max_energy=manual_deglitch_max_energy,
            manual_deglitch_margin_points=manual_deglitch_margin_points,
            save_detector_raw_overview_plot=self.plot_detector_raw_overview.get(),
            save_processed_overview_plot=self.plot_processed_overview.get(),
            save_bkgcorr_overview_plot=self.plot_bkgcorr_overview.get(),
            save_raw_overview_plot=False,
            save_norm_overview_plot=self.plot_norm_overview.get(),
            save_replicate_qc_plots=self.plot_replicate_qc.get(),
            save_drift_plot=self.save_drift_plot.get(),
            save_foil_alignment_plots=False,
            plot_energy_min=plot_min,
            plot_energy_max=plot_max,
        )

    def config_to_dict(self) -> dict:
        c = self.build_config()
        return {
            "analysis_mode": c.analysis_mode,
            "alignment_source": c.alignment_source,
            "foil_alignment_mode": c.foil_alignment_mode,
            "foil_keyword": c.foil_keyword,
            "e0": c.e0,
            "pre1": c.pre1,
            "pre2": c.pre2,
            "norm1": c.norm1,
            "norm2": c.norm2,
            "nnorm": c.nnorm,
            "align_window_min": c.align_window_min,
            "align_window_max": c.align_window_max,
            "shift_bound_min": c.shift_bound_min,
            "shift_bound_max": c.shift_bound_max,
            "alignment_quality_warn_threshold": c.alignment_quality_warn_threshold,
            "alignment_grid_points": c.alignment_grid_points,
            "fluo_multiplicative_constant": c.fluo_multiplicative_constant,
            "enable_deglitching": self.enable_deglitching.get(),
            "deglitch_mode": self.deglitch_mode.get(),
            "enable_auto_deglitch": c.enable_auto_deglitch,
            "deglitch_method": "interpolate",
            "deglitch_threshold": c.deglitch_threshold,
            "deglitch_window": c.deglitch_window,
            "deglitch_min_energy": c.deglitch_min_energy,
            "deglitch_max_energy": c.deglitch_max_energy,
            "enable_manual_deglitch_range": c.enable_manual_deglitch_range,
            "manual_deglitch_min_energy": c.manual_deglitch_min_energy,
            "manual_deglitch_max_energy": c.manual_deglitch_max_energy,
            "manual_deglitch_margin_points": c.manual_deglitch_margin_points,
            "save_detector_raw_overview_plot": getattr(c, "save_detector_raw_overview_plot", False),
            "save_processed_overview_plot": getattr(c, "save_processed_overview_plot", True),
            "save_bkgcorr_overview_plot": getattr(c, "save_bkgcorr_overview_plot", False),
            "save_norm_overview_plot": c.save_norm_overview_plot,
            "save_replicate_qc_plots": c.save_replicate_qc_plots,
            "save_drift_plot": c.save_drift_plot,
            "plot_energy_min": c.plot_energy_min,
            "plot_energy_max": c.plot_energy_max,
        }

    def apply_config_dict(self, data: dict):
        mapping = {
            "analysis_mode": self.analysis_mode,
            "alignment_source": self.alignment_source,
            "foil_alignment_mode": self.foil_mode,
            "foil_keyword": self.foil_keyword,
            "e0": self.e0,
            "pre1": self.pre1,
            "pre2": self.pre2,
            "norm1": self.norm1,
            "norm2": self.norm2,
            "nnorm": self.nnorm,
            "align_window_min": self.align_min,
            "align_window_max": self.align_max,
            "shift_bound_min": self.shift_min,
            "shift_bound_max": self.shift_max,
            "alignment_quality_warn_threshold": self.alignment_quality_warn_threshold,
            "alignment_grid_points": self.alignment_grid_points,
            "fluo_multiplicative_constant": self.fluo_factor,
            "enable_deglitching": self.enable_deglitching,
            "deglitch_mode": self.deglitch_mode,
            "enable_auto_deglitch": self.enable_auto_deglitch,
            "deglitch_threshold": self.deglitch_threshold,
            "deglitch_window": self.deglitch_window,
            "deglitch_min_energy": self.deglitch_min_energy,
            "deglitch_max_energy": self.deglitch_max_energy,
            "enable_manual_deglitch_range": self.enable_manual_deglitch_range,
            "manual_deglitch_min_energy": self.manual_deglitch_min_energy,
            "manual_deglitch_max_energy": self.manual_deglitch_max_energy,
            "manual_deglitch_margin_points": self.manual_deglitch_margin_points,
            "save_detector_raw_overview_plot": self.plot_detector_raw_overview,
            "save_processed_overview_plot": self.plot_processed_overview,
            "save_bkgcorr_overview_plot": self.plot_bkgcorr_overview,
            "save_raw_overview_plot": self.plot_processed_overview,  # legacy config support
            "save_norm_overview_plot": self.plot_norm_overview,
            "save_replicate_qc_plots": self.plot_replicate_qc,
            "save_drift_plot": self.save_drift_plot,
            "plot_energy_min": self.plot_min,
            "plot_energy_max": self.plot_max,
        }

        for key, var in mapping.items():
            if key in data:
                if isinstance(var, tk.BooleanVar):
                    var.set(bool(data[key]))
                elif data[key] is None:
                    var.set("")
                else:
                    var.set(str(data[key]))

        if "enable_deglitching" not in data or "deglitch_mode" not in data:
            auto_enabled = self.enable_auto_deglitch.get()
            manual_enabled = self.enable_manual_deglitch_range.get()
            self.enable_deglitching.set(auto_enabled or manual_enabled)
            if auto_enabled and manual_enabled:
                self.deglitch_mode.set("both")
            elif manual_enabled:
                self.deglitch_mode.set("manual")
            else:
                self.deglitch_mode.set("automatic")
        elif self.deglitch_mode.get() not in {"automatic", "manual", "both"}:
            self.deglitch_mode.set("automatic")

        self.deglitch_method.set("interpolate")
        self._sync_deglitch_backend_vars()
        self._update_deglitch_ui()
        self._update_alignment_ui()

    def save_config(self):
        try:
            data = self.config_to_dict()
        except ValueError as exc:
            messagebox.showerror("Invalid parameter", str(exc))
            return

        path = filedialog.asksaveasfilename(
            title="Save ASTRA processing config",
            defaultextension=".json",
            filetypes=[("JSON config", "*.json"), ("All files", "*.*")],
        )
        if not path:
            return

        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)

        self.log(f"Saved config: {path}")

    def load_config(self):
        path = filedialog.askopenfilename(
            title="Load ASTRA processing config",
            filetypes=[("JSON config", "*.json"), ("All files", "*.*")],
        )
        if not path:
            return

        try:
            with open(path, "r", encoding="utf-8") as f:
                data = json.load(f)
            self.apply_config_dict(data)
        except Exception as exc:
            messagebox.showerror("Could not load config", str(exc))
            return

        self.log(f"Loaded config: {path}")
        
    def clear_log(self):
        self._suppress_log_until = time.time() + 0.5

        while True:
            try:
                self._log_queue.get_nowait()
            except queue.Empty:
                break

        self.log_box.configure(state="normal")
        self.log_box.delete("1.0", tk.END)
        self.log_box.configure(state="disabled")
        self.log_box.update_idletasks()

        self.status.set("Log cleared.")

    def log(self, text):
        self._log_queue.put(str(text))

    def _drain_log_queue(self):
        try:
            while True:
                text = self._log_queue.get_nowait()

                if time.time() < self._suppress_log_until:
                    continue

                self.log_box.configure(state="normal")
                self.log_box.insert("end", text + "\n")
                self.log_box.see("end")
                self.log_box.configure(state="disabled")
        except queue.Empty:
            pass

        self.after(100, self._drain_log_queue)

    def set_running(self, running: bool):
        self._running = running

        if running:
            self.status.set("Processing…")
            self.progress.start(10)
            self.run_button.configure(state="disabled")
        else:
            self.status.set("Ready")
            self.progress.stop()
            self.run_button.configure(state="normal")

    def run_processing(self):
        if self._running:
            return

        if not self.input_dir.get().strip():
            messagebox.showerror("Missing input", "Please select an input folder.")
            return

        try:
            config = self.build_config()
        except ValueError as exc:
            messagebox.showerror("Invalid parameter", str(exc))
            return

        input_dir = self.input_dir.get().strip()
        output_dir = self.output_dir.get().strip() or None

        self.clear_log()
        self.set_running(True)
        self.log("Starting processing job…")

        def worker():
            try:
                result = process_folder(input_dir, output_dir, config=config, log=self.log)
                self.log("Finished successfully.")
                self.after(
                    0,
                    lambda: messagebox.showinfo(
                        "Finished",
                        (
                            f"Processed {result['groups_processed']} group(s).\n"
                            f"Output: {result['output_dir']}\n"
                            f"Plots: {result.get('plots_dir')}"
                        ),
                    ),
                )
            except Exception as exc:
                err = str(exc)
                self.log(f"ERROR: {err}")
                self.after(0, lambda err=err: messagebox.showerror("Processing failed", err))
            finally:
                self.after(0, lambda: self.set_running(False))

        threading.Thread(target=worker, daemon=True).start()


def main():
    app = AstraGui()
    app.mainloop()


if __name__ == "__main__":
    main()
