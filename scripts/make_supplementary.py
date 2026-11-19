"""Build the supplementary tables document."""

import os

import pandas as pd
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "manuscript", "supplementary.docx")


def setup(doc):
    s = doc.sections[0]
    s.page_width, s.page_height = Cm(21.0), Cm(29.7)
    for a in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(s, a, Cm(2.0))
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(10)
    st.paragraph_format.space_after = Pt(4)


def three_line(table):
    tblPr = table._tbl.tblPr
    borders = OxmlElement("w:tblBorders")
    for edge, val, sz in (("top", "single", 8), ("left", "nil", 0), ("bottom", "single", 8),
                          ("right", "nil", 0), ("insideH", "nil", 0), ("insideV", "nil", 0)):
        el = OxmlElement(f"w:{edge}")
        el.set(qn("w:val"), val)
        if val != "nil":
            el.set(qn("w:sz"), str(sz)); el.set(qn("w:space"), "0"); el.set(qn("w:color"), "000000")
        borders.append(el)
    tblPr.append(borders)
    for cell in table.rows[0].cells:
        tcPr = cell._tc.get_or_add_tcPr()
        tcB = OxmlElement("w:tcBorders")
        bottom = OxmlElement("w:bottom")
        bottom.set(qn("w:val"), "single"); bottom.set(qn("w:sz"), "6")
        bottom.set(qn("w:space"), "0"); bottom.set(qn("w:color"), "000000")
        tcB.append(bottom); tcPr.append(tcB)


def add_table(doc, df, caption, note=None, font=8):
    p = doc.add_paragraph()
    r = p.add_run(caption)
    r.bold = True
    r.font.size = Pt(10)
    p.paragraph_format.space_before = Pt(10)
    p.paragraph_format.keep_with_next = True

    t = doc.add_table(rows=1, cols=len(df.columns))
    t.style = "Table Grid"
    for k, col in enumerate(df.columns):
        cell = t.rows[0].cells[k]
        cell.text = ""
        run = cell.paragraphs[0].add_run(str(col))
        run.bold = True
        run.font.size = Pt(font)
    trPr = t.rows[0]._tr.get_or_add_trPr()
    el = OxmlElement("w:tblHeader")
    el.set(qn("w:val"), "true")
    trPr.append(el)
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for k, val in enumerate(row):
            cells[k].text = ""
            run = cells[k].paragraphs[0].add_run("" if pd.isna(val) else str(val))
            run.font.size = Pt(font)
    three_line(t)
    if note:
        np_ = doc.add_paragraph()
        nr = np_.add_run(note)
        nr.font.size = Pt(8)
        nr.italic = True


def main():
    doc = Document()
    setup(doc)

    h = doc.add_paragraph()
    hr = h.add_run("Supplementary material")
    hr.bold = True
    hr.font.size = Pt(14)
    sub = doc.add_paragraph()
    sr = sub.add_run(
        "Two mechanisms of lactic-acid toxicity explain pH-dependent D-lactate tolerance in "
        "Saccharomyces cerevisiae: a pH-resolved enzyme-constrained dynamic flux balance analysis"
    )
    sr.italic = True
    sr.font.size = Pt(10)

    doc.add_paragraph(
        "Table S1 lists every parameter used in the model with its value, unit and source. "
        "Table S2 lists the eight strains used for validation. Table S3 reports the propagation "
        "of a ±20% uncertainty in the two toxicity coefficients to the predicted titre."
    )

    p1 = pd.read_csv(os.path.join(ROOT, "data", "parameters.csv"))
    add_table(doc, p1, "Table S1. Model parameters.",
              "Sources: 'this work' refers to values obtained by calibration against [1]; "
              "'assumption' refers to physically plausible ranges tested in Fig. 5.")

    p2 = pd.read_csv(os.path.join(ROOT, "data", "literature_strains.csv"))
    p2 = p2[["strain", "mode", "pH_control", "titer_gL", "yield_g_g", "doi"]]
    p2.columns = ["Strain", "Mode", "pH control", "Titre (g L-1)", "Yield (g g-1)", "DOI"]
    add_table(doc, p2, "Table S2. D-lactate-producing strains used for validation.",
              "The final row reports a yield only; the corresponding titre was not reported in the source.")

    p3 = pd.read_csv(os.path.join(ROOT, "results", "param_uncertainty.csv"))
    p3["f_kh"] = p3["f_kh"].map(lambda x: f"{x:.1f}")
    p3["f_ka"] = p3["f_ka"].map(lambda x: f"{x:.1f}")
    p3["titer"] = p3["titer"].round(1)
    p3.columns = ["k_H factor", "k_A factor", "pH", "Predicted titre (g L-1)"]
    add_table(doc, p3, "Table S3. Sensitivity of the predicted titre to the two toxicity coefficients.",
              "The central row (1.0, 1.0) corresponds to the calibrated values.")

    doc.save(OUT)
    print("saved", OUT)


if __name__ == "__main__":
    main()
