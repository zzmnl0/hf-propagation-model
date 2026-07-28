"""Update the reflection-route DOCX with five-scene O/X and seeded P2P evidence."""

from __future__ import annotations

import json
import os
from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt
from docx.text.paragraph import Paragraph


ROOT = Path(__file__).resolve().parents[1]
DOCX = ROOT / "word" / "hf_3d_hf_reflection_model_route_and_progress_2026-07-26.docx"
OUT = ROOT / "output2"
ASSETS = OUT / "report_assets" / "reflection_route_2026-07-26"
ACCEPTANCE = ASSETS / "formal_complex_scene_ox_p2p_acceptance.json"
SCENES = ("baseline", "tid", "spreadf", "es", "bubble")
SCENE_NAMES = {
    "baseline": "Baseline",
    "tid": "TID",
    "spreadf": "Spread-F",
    "es": "Es",
    "bubble": "Bubble",
}


def _find_paragraph(doc, prefix):
    return next(p for p in doc.paragraphs if p.text.startswith(prefix))


def _replace_text(paragraph, text):
    paragraph.clear()
    paragraph.add_run(text)


def _add_before(doc, reference, text="", style=None, align=None):
    paragraph = doc.add_paragraph(style=style)
    if text:
        paragraph.add_run(text)
    if align is not None:
        paragraph.alignment = align
    reference._p.addprevious(paragraph._p)
    return paragraph


def _add_picture_before(doc, reference, path, width=6.35):
    paragraph = _add_before(
        doc, reference, align=WD_ALIGN_PARAGRAPH.CENTER)
    paragraph.add_run().add_picture(str(path), width=Inches(width))
    return paragraph


def _page_break_before(doc, reference):
    paragraph = _add_before(doc, reference)
    paragraph.add_run().add_break(WD_BREAK.PAGE)


def _set_cell_width(cell, width_inches):
    width = Inches(width_inches)
    cell.width = Inches(width_inches)
    tc_pr = cell._tc.get_or_add_tcPr()
    tc_w = tc_pr.find(qn("w:tcW"))
    if tc_w is None:
        tc_w = OxmlElement("w:tcW")
        tc_pr.append(tc_w)
    tc_w.set(qn("w:w"), str(width.twips))
    tc_w.set(qn("w:type"), "dxa")


def _format_table(table, widths):
    table.style = "Table Grid"
    table.autofit = False
    for row_index, row in enumerate(table.rows):
        for index, (cell, width) in enumerate(zip(row.cells, widths)):
            _set_cell_width(cell, width)
            cell.vertical_alignment = 1
            for paragraph in cell.paragraphs:
                paragraph.alignment = (
                    WD_ALIGN_PARAGRAPH.CENTER if index < 4
                    else WD_ALIGN_PARAGRAPH.LEFT)
                for run in paragraph.runs:
                    run.font.size = Pt(8.5)
                    if row_index == 0:
                        run.bold = True


def _termination_text(diagnostics):
    counts = diagnostics["termination_counts"]
    return "、".join(
        "{} {}".format(name, counts[name]) for name in sorted(counts))


def _root_text(scene, metric):
    o_count = int(metric["p2p_root_counts"]["O"])
    x_count = int(metric["p2p_root_counts"]["X"])
    if o_count + x_count == 0:
        return (
            "{}在每分支1个最近地面种子、最多5次TRF评估的受限搜索中未得到"
            "经research终值复算仍命中接收机的根。该结果只表示局部种子未收敛，"
            "不能解释为该环境不存在传播模式或不存在O/X分支。").format(
                SCENE_NAMES[scene])
    error = float(metric["p2p_max_landing_error_km"])
    return (
        "{}获得O波{}根、X波{}根，research终值复算后的最大接收机落点误差为"
        "{:.6f} km。图中轨迹是已找到的接收机命中解，不代表全角域全部模式。"
    ).format(SCENE_NAMES[scene], o_count, x_count, error)


def _replace_figure3(doc):
    caption = _find_paragraph(doc, "图3 TID")
    image_paragraph = caption._p.getprevious()
    for child in list(image_paragraph):
        image_paragraph.remove(child)
    paragraph = Paragraph(image_paragraph, caption._parent)
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    paragraph.paragraph_format.page_break_before = True
    paragraph.add_run().add_picture(
        str(ASSETS / "03_ne_relative_perturbations.png"),
        width=Inches(6.35))
    _replace_text(
        _find_paragraph(doc, "数据来源：同一正式Provider实算；各扰动统一除以"),
        "数据来源：同一正式Provider实算；四类扰动统一除以Baseline峰值。"
        "Es面板增加hEs±0.6 km显示放大窗及独立物理色标，传播Ne未缩放。")
    _replace_text(
        _find_paragraph(doc, "相对扰动图以统一Baseline峰值归一化"),
        "相对扰动图以统一Baseline峰值归一化。TID呈长尺度波动，Spread-F表现"
        "为F区随机结构，Bubble表现为F区场向耗空管；Es薄层在全高度坐标中"
        "像素占比很小，因此以110 km中心线和局部放大窗显示真实ΔNe。该放大"
        "只改变绘图视窗，Provider哈希和传播采样值逐点不变。")


def _replace_section6(doc, metrics):
    heading6 = _find_paragraph(doc, "6 复杂环境反射结果")
    heading7 = _find_paragraph(doc, "7 外部软件验证阶段简述")
    body = doc._element.body
    elements = list(body)
    start = elements.index(heading6._p)
    end = elements.index(heading7._p)
    for element in elements[start + 1:end]:
        body.remove(element)

    _add_before(
        doc, heading7,
        "在相同10 MHz、单程、球形地球和ECEF AH3D条件下，五类Provider均"
        "完成O波及O/X射线扇计算。每个波模采用33个固定发射仰角；Ray-PE、"
        "P-D和散射功率改写关闭。O波既有冻结PNG/JSON哈希保持不变，新增O/X"
        "结果另行命名，避免覆盖历史证据。")
    _add_before(
        doc, heading7,
        "五场景P2P采用质量—效率平衡配置：O/X分别选取1个最接近接收机的"
        "地面射线作为局部种子，smoke配置参与TRF优化，每根最多5次残差评估；"
        "收敛角度再以research配置前向复算，只有ground终止且落点误差小于"
        "5 km的轨迹才进入正式图表。该流程不简化Hamiltonian方程、色散投影"
        "或多尺度中心差分，但搜索范围受限，不能替代P2-6H。")

    table = doc.add_table(rows=1, cols=6)
    heading7._p.addprevious(table._tbl)
    headers = ("环境", "O根", "X根", "最大误差/km", "O终止状态", "X终止状态")
    for cell, text in zip(table.rows[0].cells, headers):
        cell.text = text
    for scene in SCENES:
        metric = metrics[scene]
        row = table.add_row().cells
        max_error = metric["p2p_max_landing_error_km"]
        values = (
            SCENE_NAMES[scene],
            str(metric["p2p_root_counts"]["O"]),
            str(metric["p2p_root_counts"]["X"]),
            "无收敛根" if max_error is None else "{:.6f}".format(max_error),
            _termination_text(metric["ray_fan"]["O"]),
            _termination_text(metric["ray_fan"]["X"]),
        )
        for cell, text in zip(row, values):
            cell.text = text
    _format_table(table, (0.7, 0.45, 0.45, 0.85, 1.9, 1.9))
    _add_before(
        doc, heading7,
        "表6-1 五类参数化环境O/X Ray Fan与局部种子P2P结果汇总",
        style="Caption Text")

    _add_before(doc, heading7, "6.1 五场景正式O波Ray Fan", style="Heading 2")
    _add_before(
        doc, heading7,
        "以下五张画布复用既有正式ECEF O波结果，统一采用33个仰角。图像背景"
        "和传播计算均来自相同场景Provider；Ray Fan表示发射角采样响应，不是"
        "接收机约束模式集合。")
    for index, scene in enumerate(SCENES):
        if index:
            _page_break_before(doc, heading7)
        _add_picture_before(
            doc, heading7, OUT / "tx_fan_unified_{}.png".format(scene))
        _add_before(
            doc, heading7,
            "图11-1{} {}环境正式ECEF AH3D O波Ray Fan".format(
                chr(ord("a") + index), SCENE_NAMES[scene]),
            style="Caption Text")
        _add_before(
            doc, heading7,
            "计算条件：10 MHz、单程、O波、33个固定仰角；P2P/P-D/Ray-PE/"
            "散射功率关闭。O波终止统计为：{}。".format(
                _termination_text(metrics[scene]["ray_fan"]["O"])))

    _page_break_before(doc, heading7)
    _add_before(
        doc, heading7, "6.2 五场景O/X局部种子P2P结果", style="Heading 2")
    _add_before(
        doc, heading7,
        "图11-2a至11-2e仅显示经research终值复算并满足接收机命中门槛的"
        "局部种子根。Spread-F与Es的空轨迹画布是显式未收敛诊断，不是无模式"
        "结论。")
    for index, scene in enumerate(SCENES):
        if index:
            _page_break_before(doc, heading7)
        _add_picture_before(
            doc, heading7, OUT / "p2p_modes_unified_{}.png".format(scene))
        _add_before(
            doc, heading7,
            "图11-2{} {}环境ECEF AH3D O/X局部种子P2P结果".format(
                chr(ord("a") + index), SCENE_NAMES[scene]),
            style="Caption Text")
        _add_before(doc, heading7, _root_text(scene, metrics[scene]))

    _page_break_before(doc, heading7)
    _add_before(
        doc, heading7, "6.3 五场景正式O/X Ray Fan", style="Heading 2")
    _add_before(
        doc, heading7,
        "每个环境在同一画布叠加O波红色实线和X波蓝色虚线，每分支各33次"
        "发射。两分支的落地、网格退出和物理域拒绝分别计数，不将拒绝轨迹"
        "静默删除。")
    for index, scene in enumerate(SCENES):
        if index:
            _page_break_before(doc, heading7)
        _add_picture_before(
            doc, heading7, OUT / "tx_fan_ox_unified_{}.png".format(scene))
        _add_before(
            doc, heading7,
            "图11-3{} {}环境正式ECEF AH3D O/X Ray Fan".format(
                chr(ord("a") + index), SCENE_NAMES[scene]),
            style="Caption Text")
        _add_before(
            doc, heading7,
            "O波：{}；X波：{}。该图用于比较磁离子分支的轨迹与终止状态，"
            "不据此推断模式完整性或绝对功率。".format(
                _termination_text(metrics[scene]["ray_fan"]["O"]),
                _termination_text(metrics[scene]["ray_fan"]["X"])))

    _page_break_before(doc, heading7)
    _add_picture_before(
        doc, heading7, ASSETS / "13_formal_ray_fan_statuses.png")
    _add_before(
        doc, heading7,
        "图11-4 五类参数化环境冻结O波Ray Fan终止状态统计",
        style="Caption Text")
    _add_before(
        doc, heading7,
        "图11-4保留原O波冻结统计，用于与拆分后的图11-1a至11-1e对应。"
        "新增O/X逐分支统计见表6-1及机器可读验收JSON。")


def _update_other_sections(doc, metrics):
    _replace_text(
        _find_paragraph(doc, "统一场景入口现已闭合地面至60 km"),
        "统一场景入口现已闭合地面至60 km真空域、60—600 km数据域和严格"
        "水平边界。五类环境均已生成冻结O波Ray Fan、新增O/X Ray Fan和受限"
        "局部种子P2P图。传播计算与背景可视化来自同一Provider；P-D、Ray-PE"
        "和散射功率改写保持关闭。")
    p2_limited = _find_paragraph(doc, "当前仍记录52个未进一步正式求解")
    new_paragraph = OxmlElement("w:p")
    p2_limited._p.addnext(new_paragraph)
    paragraph = Paragraph(new_paragraph, p2_limited._parent)
    paragraph.add_run(
        "本报告新增的五场景P2P属于不同语义：仅以每分支1个Ray Fan近邻种子"
        "执行受限局部求根，并以research前向复算确认终值。其结果用于展示已"
        "找到的接收机命中路径，不能覆盖P2-6L参考算例，也不能替代延期的"
        "P2-6H复杂场全角域完整性审计。")

    for paragraph in doc.paragraphs:
        if paragraph.text.startswith("保持已闭合的地面真空域"):
            _replace_text(
                paragraph,
                "保持已闭合的地面真空域、PCHIP垂直插值、五场景冻结O波结果及"
                "本次新增O/X Ray Fan和局部种子P2P证据；后续若增加种子、扩大"
                "角域或提高TRF上限，必须形成新版本，不静默覆盖当前有限结果。")
        elif paragraph.text.startswith("阶段工作已形成复杂电离层三维介质输入"):
            _replace_text(
                paragraph,
                "阶段工作已形成复杂电离层三维介质输入、地面真空域、AH3D"
                "磁离子反射、高低射线参考识别和内部精度证据链；五类参数化环境"
                "已完成冻结O波Ray Fan、O/X Ray Fan及受限局部种子P2P可视化。"
                "Baseline和Bubble获得O/X命中根，TID获得O根；Spread-F和Es在"
                "当前受限配置下未收敛，不能据此判定无传播模式。")

    scene_table = doc.tables[5]
    current = {
        "Baseline": "O及O/X Ray Fan；局部P2P获得O/X根",
        "TID": "O及O/X Ray Fan；局部P2P获得O根",
        "Spread-F": "O及O/X Ray Fan；局部种子未收敛",
        "Es": "O及O/X Ray Fan；局部种子未收敛",
        "Bubble": "O及O/X Ray Fan；局部P2P获得O/X根",
    }
    for row in scene_table.rows[1:]:
        name = row.cells[0].text
        if name in current:
            row.cells[1].text = current[name]
            row.cells[2].text = (
                "比较已找到根的落点、反射高度、路径量与群时延；"
                "未收敛场景需扩大种子/角域后复验")
            row.cells[3].text = "局部种子结果不等于复杂场全模式完整性"

    evidence = doc.tables[7]
    row = evidence.add_row().cells
    values = (
        "复杂场O/X与局部P2P",
        "output2/tx_fan_ox_unified_*、p2p_modes_unified_*；"
        "output2/report_assets/reflection_route_2026-07-26",
        "PNG、CSV/JSON、日志、SHA-256 manifest",
        "受限局部种子有限结果；P2-6H仍未验证",
    )
    for cell, text in zip(row, values):
        cell.text = text


def run():
    with ACCEPTANCE.open("r", encoding="utf-8") as stream:
        metrics = json.load(stream)["metrics"]
    doc = Document(str(DOCX))
    _replace_figure3(doc)
    _update_other_sections(doc, metrics)
    _replace_section6(doc, metrics)
    temporary = DOCX.with_suffix(".updated.docx")
    doc.save(str(temporary))
    os.replace(temporary, DOCX)
    print(DOCX)


if __name__ == "__main__":
    run()
