## Project pipeline diagram

```mermaid
flowchart TD

    A[Raw microscopy dataset<br/>TIFF stack<br/>3 channels × 121 time points] --> B[01_organize_images.py<br/>Split stack by channel and frame]

    B --> C1[Phase contrast folder<br/>PC frames 001-121]
    B --> C2[ZapA folder<br/>ZapA frames 001-121]
    B --> C3[PBP2B folder<br/>PBP2B frames 001-121]

    C1 --> D1[02_run_pretrained_cellpose.py<br/>Pretrained Cellpose v3 on phase contrast]
    C3 --> D2[03_preprocess_signal_channels.py<br/>Noise reduction on PBP2B]

    D2 --> D3[Test Cellpose on denoised PBP2B]
    D1 --> E[Segmentation quality check<br/>Masks, outlines, screenshots]
    D3 --> E

    E --> F[04_prepare_training_data.py<br/>Select training frames<br/>Control and sample datasets]

    F --> G[Manual correction<br/>ImageJ or napari<br/>Correct Cellpose masks]

    G --> H[05_train_cellpose_model.py<br/>Train custom Cellpose v3 model]

    H --> I[Validate trained model<br/>Compare predicted masks with corrected masks]

    I --> J{Are masks good enough?}

    J -- No --> F
    J -- Yes --> K[06_segment_full_dataset.py<br/>Segment all frames with trained model]

    K --> L[Final segmentation outputs<br/>Masks, outlines, GIFs]

    L --> M[07_track_cells.py<br/>Track cells over time<br/>LapTrack or TrackMate]

    M --> N[08_measure_cells_and_signals.py<br/>Measure shape and fluorescence]

    N --> O1[Cell length over time]
    N --> O2[ZapA mid-cell intensity over time]
    N --> O3[PBP2B mid-cell intensity over time]
    N --> O4[Division timing]

    O1 --> P[09_plot_results.py<br/>Create plots]
    O2 --> P
    O3 --> P
    O4 --> P

    P --> Q[Final outputs<br/>CSV files, plots, GIFs, trained models]

    Q --> R[README.md report<br/>Results, validation, interpretation]
    Q --> S[Presentation<br/>Background, approach, issues, results, outlook]
