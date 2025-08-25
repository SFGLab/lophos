from pathlib import Path

import pandas as pd

def write_summary(path: Path, peaks: pd.DataFrame, loops: pd.DataFrame) -> None:
    s = []
    s.append(("peaks_total", len(peaks)))
    s.append(("peaks_maternal", int((peaks["bias_call"] == "Maternal").sum())))
    s.append(("peaks_paternal", int((peaks["bias_call"] == "Paternal").sum())))
    s.append(("peaks_balanced", int((peaks["bias_call"] == "Balanced").sum())))
    s.append(("peaks_undetermined", int((peaks["bias_call"] == "Undetermined").sum())))
    s.append(("loops_total", len(loops)))
    s.append(("loops_maternal", int((loops["bias_call"] == "Maternal").sum())))
    s.append(("loops_paternal", int((loops["bias_call"] == "Paternal").sum())))
    s.append(("loops_balanced", int((loops["bias_call"] == "Balanced").sum())))
    s.append(("loops_undetermined", int((loops["bias_call"] == "Undetermined").sum())))
    pd.DataFrame(s, columns=["metric", "value"]).to_csv(path, sep="\t", index=False)
