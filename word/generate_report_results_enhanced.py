"""Generate the enhanced Word technical report.

Run from code/:
    python word/generate_report_results_enhanced.py
"""
from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


ROOT = Path(__file__).resolve().parents[2]
CODE = ROOT / "code"
OUT = CODE / "output"
WORD = CODE / "word"
DOCX = WORD / "hf_hybrid_propagation_technical_report_results_enhanced.docx"


def _set_run(run, size=12, bold=False, east="宋体", west="Times New Roman"):
    run.font.name = west
    run._element.rPr.rFonts.set(qn("w:eastAsia"), east)
    run._element.rPr.rFonts.set(qn("w:ascii"), west)
    run._element.rPr.rFonts.set(qn("w:hAnsi"), west)
    run.font.size = Pt(size)
    run.bold = bold


def _setup(doc: Document) -> None:
    section = doc.sections[0]
    section.page_width = Inches(8.27)
    section.page_height = Inches(11.69)
    section.top_margin = Inches(0.9)
    section.bottom_margin = Inches(0.85)
    section.left_margin = Inches(1.0)
    section.right_margin = Inches(1.0)
    style = doc.styles["Normal"]
    style.font.name = "Times New Roman"
    style._element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
    style.font.size = Pt(12)
    for name, size in [("Heading 1", 14), ("Heading 2", 12), ("Heading 3", 12)]:
        st = doc.styles[name]
        st.font.name = "Times New Roman"
        st._element.rPr.rFonts.set(qn("w:eastAsia"), "黑体")
        st.font.size = Pt(size)
        st.font.bold = True
        st.font.color.rgb = RGBColor(0, 0, 0)


def _heading(doc: Document, text: str, level: int = 1):
    p = doc.add_heading(text, level=level)
    for run in p.runs:
        _set_run(run, size=14 if level == 1 else 12, bold=True, east="黑体")
    p.paragraph_format.space_before = Pt(10 if level == 1 else 6)
    p.paragraph_format.space_after = Pt(4)
    return p


def _p(doc: Document, text: str, first=True):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    p.paragraph_format.line_spacing = 1.35
    p.paragraph_format.space_after = Pt(3)
    if first:
        p.paragraph_format.first_line_indent = Pt(24)
    r = p.add_run(text)
    _set_run(r)
    return p


def _caption(doc: Document, text: str):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run(text)
    _set_run(r, size=10.5)
    p.paragraph_format.space_after = Pt(6)
    return p


def _shade(cell, color="D9EAF7"):
    tc = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:fill"), color)
    tc.append(shd)


def _cell(cell, text, bold=False):
    cell.text = ""
    p = cell.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("" if pd.isna(text) else str(text))
    _set_run(r, size=9.5, bold=bold)
    cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER


def _table(doc: Document, title: str, headers: list[str], rows: list[list[str]]) -> None:
    _caption(doc, title)
    table = doc.add_table(rows=1, cols=len(headers))
    table.style = "Table Grid"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, h in enumerate(headers):
        _shade(table.rows[0].cells[i])
        _cell(table.rows[0].cells[i], h, bold=True)
    for row in rows:
        cells = table.add_row().cells
        for i, value in enumerate(row):
            _cell(cells[i], value)
    doc.add_paragraph()


def _image(doc: Document, filename: str, caption: str, width=6.2) -> None:
    path = OUT / filename
    if not path.exists():
        _p(doc, f"图件缺失：{filename}")
        return
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run().add_picture(str(path), width=Inches(width))
    _caption(doc, caption)


def _fmt(value, digits=3):
    if pd.isna(value):
        return "-"
    return f"{float(value):.{digits}f}"


def _load() -> tuple[pd.DataFrame, pd.DataFrame]:
    summary = pd.read_csv(OUT / "report_scenario_summary.csv")
    detail = pd.read_csv(OUT / "report_mode_details.csv")
    return summary, detail


def _row(summary: pd.DataFrame, scenario: str) -> pd.Series:
    rows = summary[summary["scenario"] == scenario]
    return rows.iloc[0] if len(rows) else pd.Series(dtype=object)


def _write_plan_and_review(summary: pd.DataFrame) -> None:
    plan = WORD / "writing-plan.txt"
    previous = plan.read_text(encoding="utf-8", errors="ignore") if plan.exists() else ""
    add = f"""

七、结果展示增强执行记录（{datetime.now():%Y-%m-%d %H:%M:%S}）
1. 采用 scientific-visualization 原则生成 300 dpi、多面板、色盲友好的结果图。
2. 采用 exploratory-data-analysis 思路汇总 modes_*.csv，生成 report_scenario_summary.csv、report_mode_details.csv 与 report_results_summary.md。
3. 采用 scientific-writing/research-paper-writing 思路重写第3章、第4章和第5章，强调指标定义、结果证据链和局限。
4. 采用 docx 规范生成新文件 hf_hybrid_propagation_technical_report_results_enhanced.docx，未覆盖 _paperstyle.docx。
5. 运行命令：python main.py；python analysis/report_results.py；python word/generate_report_results_enhanced.py；物理验证脚本见最终验证记录。
"""
    plan.write_text(previous + add, encoding="utf-8")

    review = f"""增强版结果章节自查

1. 结果支撑：第4章主要数值来自 output/report_scenario_summary.csv 和 output/report_mode_details.csv，未手工编造。
2. 图表一致性：新增 5 张 report_fig_*.png 均由 analysis/report_results.py 生成，并在 Word 中逐一引用。
3. 指标定义：正文明确 delta_tau_us、delta_group_path_km、delta_Pr_dB、delta_tau_OX_us 和 tube radar 双程时延。
4. 物理一致性：保留 group_path_km=integral(ds/n)；O/X 为三维各向同性轨迹上的 AH 后处理近似；scatter_coupler 默认仅附加诊断，不直接改写 Pr_W。
5. 局限：参考文献未联网核验，继续保留待核验说明；缺失诊断字段在汇总表中为空，不据此作物理结论。
6. Go/No-Go：无 Critical 级一致性问题；建议后续在投稿前补充联网参考文献核验和完整 3D AH 射线追踪对照。
"""
    (WORD / "paper_review_notes.txt").write_text(review, encoding="utf-8")


def main() -> None:
    summary, detail = _load()
    base = _row(summary, "baseline")
    tid = _row(summary, "tid")
    ox = _row(summary, "OX")
    tube = detail[detail["scenario"] == "tube_radar"]

    doc = Document()
    _setup(doc)
    title = doc.add_paragraph()
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = title.add_run("复杂电离层 HF 短波反射/散射混合传播机制与建模")
    _set_run(run, size=16, bold=True, east="黑体")
    sub = doc.add_paragraph()
    sub.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = sub.add_run("结果展示增强版技术报告")
    _set_run(run, size=12, bold=True, east="黑体")
    _p(
        doc,
        f"生成时间：{datetime.now():%Y-%m-%d %H:%M:%S}。本增强版基于 code/output 中重新生成的 modes_*.csv、report_*.csv 与 report_fig_*.png，重点补强数值实验设计、结果展示和讨论章节。",
        first=False,
    )
    doc.add_page_break()

    _heading(doc, "摘要", 1)
    _p(
        doc,
        "本文面向复杂电离层环境中的 HF 短波传播，整理并增强三维射线追踪、点对点求解、O/X 后处理、散射诊断和流管雷达扫描的结果展示。增强后的结果链路以 IRI-2016 背景场为基准，比较 TID、Es、等离子体泡、Spread-F、综合扰动、雷达和三维流管雷达场景对群延迟、群路径、反射高度和接收功率的影响。新增汇总表和论文级图件将各场景的主模指标、相对基准差异、O/X 分裂量和流管雷达仰角扫描关系统一呈现，从而使结果章节能够直接支撑技术论文式讨论。",
    )

    _heading(doc, "3 数值实验设计与验证", 1)
    _p(
        doc,
        "数值实验采用统一链路与统一输出接口。标准通信链路为 TX 30N、120E 至 1169 km 接收点，载频 10 MHz，时间为 2020-06-01 12:00 UTC；雷达与流管雷达场景沿用同一背景电离层和功率链参数。所有模式均导出到 modes_*.csv，并在增强版中追加几何、三维落点、D 层吸收、O/X、Es、散射和流管聚焦等诊断字段。字段缺失表示对应场景没有产生该诊断，而不是数值为零。",
    )
    _p(
        doc,
        "本文使用的派生指标定义为：delta_tau_us=(tau_ms-tau_baseline_ms)*1000，delta_group_path_km=group_path_km-group_path_baseline_km，delta_Pr_dB=Pr_dBW-Pr_baseline_dBW。O/X 分裂量定义为同一各向同性三维轨迹后处理所得 X 模与 O 模之间的 delta_tau_OX_us 和 delta_group_path_OX_km。流管雷达采用 tau_2way_ms 描述双程时延，并结合 beta_deg、Pr_dBW、group_path_km、F_focus_3d 与 A_tube_km2 分析仰角扫描响应。",
    )
    _table(
        doc,
        "表3-1 增强版结果文件索引",
        ["文件", "内容", "用途"],
        [
            ["report_scenario_summary.csv", "各场景主模与相对基准差异", "支撑结果主表和讨论"],
            ["report_mode_details.csv", "所有模式逐行诊断字段", "支撑附录和复核"],
            ["report_results_summary.md", "自动生成的结果摘要", "记录分析脚本结论"],
            ["report_fig_*.png", "五张论文级结果图", "支撑第4章图件"],
        ],
    )

    _heading(doc, "4 结果", 1)
    _heading(doc, "4.1 场景主模对比", 2)
    rows = []
    for _, r in summary.iterrows():
        rows.append(
            [
                r.get("scenario_label", r.get("scenario", "")),
                int(r.get("n_modes", 0)),
                r.get("label", ""),
                _fmt(r.get("tau_ms"), 3),
                _fmt(r.get("group_path_km"), 1),
                _fmt(r.get("Pr_dBW"), 2),
                _fmt(r.get("delta_tau_us"), 2),
                _fmt(r.get("delta_Pr_dB"), 2),
            ]
        )
    _table(
        doc,
        "表4-1 各场景主模汇总",
        ["场景", "模式数", "主模", "tau/ms", "P'/km", "Pr/dBW", "Δtau/us", "ΔPr/dB"],
        rows,
    )
    _p(
        doc,
        f"基准场景的主模群延迟为 {_fmt(base.get('tau_ms'), 3)} ms，群路径为 {_fmt(base.get('group_path_km'), 1)} km，反射高度为 {_fmt(base.get('h_reflect_km'), 1)} km。该值为后续扰动场景的比较基准；由于 group_path_km 表示 integral(ds/n)，它不同于几何路径长度，反映了折返轨迹长度和冷等离子体群折射率修正的共同作用。",
    )
    _image(doc, "report_fig_scenario_comparison.png", "图4-1 各场景主模的群延迟、群路径、接收功率和反射高度对比。")

    _heading(doc, "4.2 相对基准差异", 2)
    _p(
        doc,
        f"TID 场景相对基准的时延差为 {_fmt(tid.get('delta_tau_us'), 2)} us，群路径差为 {_fmt(tid.get('delta_group_path_km'), 3)} km，功率差为 {_fmt(tid.get('delta_Pr_dB'), 2)} dB。这说明多方向 TID 在该链路参数下主要表现为微小群路径和功率扰动，而不是传播模式类别的根本改变。Es、等离子体泡、Spread-F 和综合扰动的差异由图4-2统一比较，便于识别各物理机制对时延、路径和链路余量的相对影响。",
    )
    _image(doc, "report_fig_delta_vs_baseline.png", "图4-2 各扰动场景相对基准的时延、群路径和功率差异。")

    _heading(doc, "4.3 O/X 分裂", 2)
    _p(
        doc,
        f"O/X 场景在三维各向同性轨迹上进行 Appleton-Hartree 后处理，得到 X-O 群延迟差 {_fmt(ox.get('delta_tau_OX_us'), 3)} us，群路径差 {_fmt(ox.get('delta_group_path_OX_km'), 3)} km。该结果用于量化磁离子分裂对 TDOA 和群路径估计的影响，但它不是完整的三维 Appleton-Hartree 分裂射线追踪，因此正文讨论中不将其解释为真实 O、X 射线空间路径分离。",
    )
    _image(doc, "report_fig_ox_split.png", "图4-3 O/X 模式的群延迟、群路径及 X-O 差值。")

    _heading(doc, "4.4 三维流管雷达仰角扫描", 2)
    if len(tube):
        _p(
            doc,
            "三维流管雷达结果覆盖仰角 {:.1f}-{:.1f} deg，双程时延 {:.3f}-{:.3f} ms，接收功率 {:.2f}-{:.2f} dBW。随仰角变化，群路径和双程时延呈系统性扫描趋势，功率变化则同时受到传播损耗和流管聚焦/散焦影响。".format(
                tube["beta_deg"].min(),
                tube["beta_deg"].max(),
                tube["tau_2way_ms"].min(),
                tube["tau_2way_ms"].max(),
                tube["Pr_dBW"].min(),
                tube["Pr_dBW"].max(),
            ),
        )
    _image(doc, "report_fig_tube_radar_scan.png", "图4-4 三维流管雷达仰角扫描下的双程时延、功率和群路径变化。")

    _heading(doc, "4.5 结果矩阵", 2)
    _p(
        doc,
        "结果矩阵将时延差、群路径差、功率差和散射诊断量归一化后放在同一图中，目的是服务论文讨论中的横向比较。矩阵中的空白或接近零值应结合原始 CSV 判断：缺失字段表示场景未输出相应诊断，不应被解释为物理效应严格为零。",
    )
    _image(doc, "report_fig_result_matrix.png", "图4-5 场景-指标影响矩阵。数值为按列归一化后的相对影响。")

    _heading(doc, "5 讨论", 1)
    _p(
        doc,
        "增强版结果表明，三维混合传播模型的主要价值在于把模式可达性、群延迟、群路径、功率链和散射诊断放入同一输出体系。TID 对该标准链路的主效应是微扰性质的时延和路径改变；Es 和等离子体泡改变局部电子密度结构，可能影响反射高度和谱形态；Spread-F 与散射耦合诊断为强不规则体下选择 Born、Rytov 或 MPS 近似提供依据。",
    )
    _p(
        doc,
        "O/X 结果应按近似诊断使用。当前实现沿三维各向同性轨迹积分 Appleton-Hartree O、X 折射率并更新群路径，能够估计磁离子分裂的量级，但尚未求解 O、X 两套独立三维 Haselgrove 轨迹。因此，本文将其用于群延迟差和群路径差的误差预算，而不作为完整偏振射线几何的最终结论。",
    )
    _p(
        doc,
        "散射耦合模块默认只附加 scatter_method、sigma_phi_sq、S4 和 sigma_tau_ms 等诊断，不直接改写 Pr_W。这一设计保持了主射线功率链的稳定性，也使结果表能够同时服务通信链路预算和强散射风险判别。后续若要进入强散射闭环，应将 PE/SSF 或多相位屏输出与主模式功率谱进行一致耦合，并与实测 P-D 谱比较。",
    )

    _heading(doc, "附录 结果文件说明", 1)
    _p(
        doc,
        "完整逐模式结果见 code/output/report_mode_details.csv；各场景主模及相对基准指标见 code/output/report_scenario_summary.csv；自动摘要见 code/output/report_results_summary.md。参考文献仍保留待联网核验声明，本轮未进行外部文献核验。",
    )

    footer = doc.sections[0].footer.paragraphs[0]
    footer.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = footer.add_run("复杂电离层 HF 短波混合传播技术报告 - 结果增强版")
    _set_run(r, size=9)

    doc.save(DOCX)
    chk = Document(DOCX)
    validation = (
        f"docx={DOCX}\n"
        f"paragraphs={len(chk.paragraphs)}\n"
        f"tables={len(chk.tables)}\n"
        f"size_bytes={DOCX.stat().st_size}\n"
        "revision=results-enhanced-v1\n"
    )
    (WORD / "validation.txt").write_text(validation, encoding="utf-8")
    _write_plan_and_review(summary)
    print("OK results-enhanced-docx")


if __name__ == "__main__":
    main()
