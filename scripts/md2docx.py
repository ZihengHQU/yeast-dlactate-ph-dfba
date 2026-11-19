"""Markdown -> DOCX converter for the manuscript (journal-submission layout)."""

import json
import os
import re
from xml.sax.saxutils import escape

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFS = json.load(open(os.path.join(ROOT, "data", "references.json"), encoding="utf-8"))
CITE = re.compile(r"\[(\d+(?:\s*[-,]\s*\d+)*)\]")
USE_FIELDS = True
ZH = False

INLINE = re.compile(r"(\*\*[^*]+\*\*|\*[^*]+\*|_\{[^}]*\}|\^\{[^}]*\}|`[^`]+`)")


def add_runs(par, text, fields=True):
    """Render **bold**, *italic*, _{sub}, ^{sup} and EndNote citation fields."""
    if fields and USE_FIELDS and CITE.search(text):
        pos = 0
        for m in CITE.finditer(text):
            if m.start() > pos:
                add_runs(par, text[pos:m.start()], fields=False)
            cite_field(par, m.group(1), m.group(0))
            pos = m.end()
        if pos < len(text):
            add_runs(par, text[pos:], fields=False)
        return
    for tok in INLINE.split(text):
        if not tok:
            continue
        if tok.startswith("**") and tok.endswith("**") and len(tok) > 4:
            r = par.add_run(tok[2:-2])
            r.bold = True
        elif tok.startswith("_{") and tok.endswith("}"):
            r = par.add_run(tok[2:-1])
            r.font.subscript = True
        elif tok.startswith("^{") and tok.endswith("}"):
            r = par.add_run(tok[2:-1])
            r.font.superscript = True
        elif tok.startswith("`") and tok.endswith("`") and len(tok) > 2:
            r = par.add_run(tok[1:-1])
            r.font.name = "Courier New"
            r.font.size = Pt(9)
        elif tok.startswith("*") and tok.endswith("*") and len(tok) > 2:
            r = par.add_run(tok[1:-1])
            r.italic = True
        else:
            par.add_run(tok)



def _expand(nums):
    out = []
    for part in nums.split(","):
        part = part.strip()
        if "-" in part:
            a, b = part.split("-")
            out.extend(range(int(a), int(b) + 1))
        elif part:
            out.append(int(part))
    return out


def cite_field(par, nums, display):
    """Insert an EndNote (ADDIN EN.CITE) field displaying the bracket text."""
    cites = ""
    for n in _expand(nums):
        r = REFS.get(str(n))
        if not r:
            continue
        title = r["raw"].split(". ")[1] if ". " in r["raw"] else ""
        cites += ("<Cite><Author>%s</Author><Year>%s</Year><RecNum>%d</RecNum>"
                  "<record><rec-number>%d</rec-number><ref-type name=\"Journal Article\">17</ref-type>"
                  "<contributors><authors><author>%s</author></authors></contributors>"
                  "<titles><title>%s</title></titles><dates><year>%s</year></dates>"
                  "<volume>%s</volume><pages>%s</pages>"
                  "<electronic-resource-num>%s</electronic-resource-num></record></Cite>") % (
            escape(r["authors"].split(",")[0]), escape(r["year"]), n, n,
            escape(r["authors"]), escape(title), escape(r["year"]),
            escape(r["volume"]), escape(r["pages"]), escape(r["doi"]))
    instr = " ADDIN EN.CITE <EndNote>%s</EndNote> " % cites
    fld = OxmlElement("w:fldSimple")
    fld.set(qn("w:instr"), instr)
    run = OxmlElement("w:r")
    t = OxmlElement("w:t")
    t.text = display
    t.set(qn("xml:space"), "preserve")
    run.append(t)
    fld.append(run)
    par._p.append(fld)


def reflist_open(par):
    run = OxmlElement("w:r")
    b = OxmlElement("w:fldChar"); b.set(qn("w:fldCharType"), "begin")
    it = OxmlElement("w:instrText"); it.set(qn("xml:space"), "preserve")
    it.text = " ADDIN EN.REFLIST "
    sp = OxmlElement("w:fldChar"); sp.set(qn("w:fldCharType"), "separate")
    for el in (b, it, sp):
        run.append(el)
    pPr = par._p.pPr
    if pPr is not None:
        pPr.addnext(run)
    else:
        par._p.insert(0, run)


def reflist_close(par):
    r = par.add_run()._r
    e = OxmlElement("w:fldChar"); e.set(qn("w:fldCharType"), "end")
    r.append(e)


def page_setup(doc):
    section = doc.sections[0]
    section.page_width = Cm(21.0)
    section.page_height = Cm(29.7)
    for attr in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(section, attr, Cm(2.0))

    ln = OxmlElement("w:lnNumType")
    ln.set(qn("w:countBy"), "1")
    ln.set(qn("w:restart"), "continuous")
    ln.set(qn("w:distance"), "360")
    section._sectPr.append(ln)

    par = section.footer.paragraphs[0]
    par.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = par.add_run()
    begin = OxmlElement("w:fldChar")
    begin.set(qn("w:fldCharType"), "begin")
    instr = OxmlElement("w:instrText")
    instr.set(qn("xml:space"), "preserve")
    instr.text = " PAGE "
    end = OxmlElement("w:fldChar")
    end.set(qn("w:fldCharType"), "end")
    for el in (begin, instr, end):
        run._r.append(el)



def three_line(table):
    """Booktabs-style rules: top, header bottom, table bottom; no vertical lines."""
    tblPr = table._tbl.tblPr
    borders = OxmlElement("w:tblBorders")
    for edge, val, sz in (("top", "single", 8), ("left", "nil", 0), ("bottom", "single", 8),
                          ("right", "nil", 0), ("insideH", "nil", 0), ("insideV", "nil", 0)):
        el = OxmlElement(f"w:{edge}")
        el.set(qn("w:val"), val)
        if val != "nil":
            el.set(qn("w:sz"), str(sz))
            el.set(qn("w:space"), "0")
            el.set(qn("w:color"), "000000")
        borders.append(el)
    tblPr.append(borders)

    for cell in table.rows[0].cells:
        tcPr = cell._tc.get_or_add_tcPr()
        tcBorders = OxmlElement("w:tcBorders")
        bottom = OxmlElement("w:bottom")
        bottom.set(qn("w:val"), "single")
        bottom.set(qn("w:sz"), "6")
        bottom.set(qn("w:space"), "0")
        bottom.set(qn("w:color"), "000000")
        tcBorders.append(bottom)
        tcPr.append(tcBorders)


def repeat_header(table):
    trPr = table.rows[0]._tr.get_or_add_trPr()
    el = OxmlElement("w:tblHeader")
    el.set(qn("w:val"), "true")
    trPr.append(el)


def is_meta_line(text):
    return text.startswith(("**Authors", "**Corresponding", "**Target", "**Keywords"))


def convert(md_path, docx_path, zh=False, use_fields=True):
    global ZH, USE_FIELDS
    ZH, USE_FIELDS = zh, use_fields
    lines = open(md_path, encoding="utf-8").read().split("\n")
    doc = Document()
    page_setup(doc)

    normal = doc.styles["Normal"]
    normal.font.name = "Times New Roman"
    if ZH:
        rpr = normal.element.get_or_add_rPr()
        rfonts = rpr.find(qn("w:rFonts"))
        if rfonts is None:
            rfonts = OxmlElement("w:rFonts"); rpr.append(rfonts)
        rfonts.set(qn("w:eastAsia"), "SimSun")
    normal.font.size = Pt(11)
    normal.paragraph_format.space_after = Pt(6)
    normal.paragraph_format.line_spacing = 1.5

    i, n = 0, len(lines)
    while i < n:
        ln = lines[i].rstrip()
        if not ln.strip():
            i += 1
            continue

        if ln.startswith("|") and i + 1 < n and set(lines[i + 1].replace("|", "").strip()) <= set("-: "):
            header = [c.strip() for c in ln.strip().strip("|").split("|")]
            i += 2
            rows = []
            while i < n and lines[i].startswith("|"):
                rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
                i += 1
            if doc.paragraphs:
                doc.paragraphs[-1].paragraph_format.keep_with_next = True
            t = doc.add_table(rows=1, cols=len(header))
            t.style = "Table Grid"
            t.alignment = WD_TABLE_ALIGNMENT.CENTER
            t.autofit = True
            for k, h in enumerate(header):
                cell = t.rows[0].cells[k]
                cell.text = ""
                par = cell.paragraphs[0]
                add_runs(par, h)
                for r in par.runs:
                    r.font.size = Pt(9)
                    r.bold = True
            widths = {6: [3.8, 5.0, 2.0, 2.4, 2.0, 1.8], 4: [5.5, 3.5, 4.0, 4.0]}.get(len(header))
            if widths:
                t.autofit = False
                for k, wd in enumerate(widths):
                    t.columns[k].width = Cm(wd)
            repeat_header(t)
            for row in rows:
                cells = t.add_row().cells
                for k, c in enumerate(row[:len(header)]):
                    cells[k].text = ""
                    par = cells[k].paragraphs[0]
                    add_runs(par, c)
                    for r in par.runs:
                        r.font.size = Pt(9)
            three_line(t)
            for row in t.rows[:-1]:
                for cell in row.cells:
                    for par in cell.paragraphs:
                        par.paragraph_format.keep_with_next = True
            doc.add_paragraph()
            continue

        if ln.startswith("!["):
            m = re.match(r"!\[[^\]]*\]\(([^)]+)\)", ln.strip())
            if m:
                img = os.path.join(ROOT, m.group(1))
                tif = img[:-4] + ".tif"
                if os.path.exists(tif):
                    img = tif
                if os.path.exists(img):
                    width = Cm(10.0) if "fig3" in img else (Cm(13.0) if "fig5" in img else Cm(15.0))
                    doc.add_picture(img, width=width)
                    doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
                    doc.paragraphs[-1].paragraph_format.keep_with_next = True
                    doc.paragraphs[-1].paragraph_format.space_before = Pt(10)
            i += 1
            continue

        if ln.startswith("#"):
            lvl = len(ln) - len(ln.lstrip("#"))
            h = doc.add_heading(level=min(lvl, 4))
            add_runs(h, ln.lstrip("#").strip())
            for r in h.runs:
                r.font.color.rgb = RGBColor(0, 0, 0)
            i += 1
            continue

        if ln.strip() == "---":
            i += 1
            continue

        if ln.startswith("  "):
            par = doc.add_paragraph()
            par.paragraph_format.left_indent = Cm(1)
            par.paragraph_format.space_after = Pt(2)
            par.paragraph_format.line_spacing = 1.0
            add_runs(par, ln.strip())
            for r in par.runs:
                r.font.name = "Courier New"
                r.font.size = Pt(9)
            i += 1
            continue

        par = doc.add_paragraph()
        add_runs(par, ln.strip())
        par.alignment = WD_ALIGN_PARAGRAPH.LEFT if is_meta_line(ln) else WD_ALIGN_PARAGRAPH.JUSTIFY
        i += 1

    refs = [p for p in doc.paragraphs if re.match(r"^\d{1,2}\.\s", p.text)]
    if refs:
        reflist_open(refs[0])
        reflist_close(refs[-1])

    doc.save(docx_path)
    print("saved", docx_path)


if __name__ == "__main__":
    import sys
    flags = [a for a in sys.argv[1:] if a.startswith("--")]
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    src = args[0] if args else os.path.join(ROOT, "manuscript", "manuscript.md")
    dst = args[1] if len(args) > 1 else os.path.join(ROOT, "manuscript", "manuscript.docx")
    convert(src, dst, zh="--zh" in flags, use_fields="--no-fields" not in flags)
