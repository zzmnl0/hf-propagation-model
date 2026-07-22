from pathlib import Path

from docx import Document
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.shared import Pt


WORD_DIR = Path(__file__).resolve().parent
OUTPUT = WORD_DIR / "hf_3d_reflection_scattering_weekly_progress_2026-07-17.docx"


LAST_WEEK = [
    (
        "统一模型与场景入口",
        "建立统一 NeProvider3D 和 main2.py pipeline，使 baseline、TID、Spread-F、Es、bubble 等参数化三维电离层场调用同一传播程序；main.py 保留单场景验证，输出语义区分单程/双程、O/X、反射/散射状态。",
    ),
    (
        "三维磁离子反射主链",
        "完成 AH3D 有限差分 Hamiltonian first version 及 V2 入口替换，支持 O/X 分支、三维射线扇、P2P、群路径和群时延；ray fan 默认仅算 O 波，可开关启用并区分 O/X，P2P 保留 O/X。",
    ),
    (
        "反射/散射与观测输出",
        "形成 Ray-PE local_2p5d、路径相关 MPS 统计核、复折射率吸收及通信/雷达 P-D 输出链；完成复杂电离层背景绘图规范和主要输出语义标注。",
    ),
    (
        "阶段边界",
        "明确 local_2p5d 不是完整 3D PE/全波，MPS 不是完整三维散射场，AH3D 仍需加强截止区、步长与跨频连续性；绝对功率尚无实测校准。",
    ),
]


THIS_WEEK = [
    (
        "A3Rnew 路线调整",
        "将 P2P 升级路线调整为精度优先：生产求解采用 DOP853、自适应积分、色散面投影和多尺度中心差分 Jacobian；变分 Jacobian 降为实验诊断，不再作为 P2/A4 准入条件。",
        "已完成",
    ),
    (
        "P0C 高低射线参考",
        "基于 Coleman 直接变分方法及 Nosikov 广义力方法整理并实现 Chapman F2 高/低射线参考，保留高射线极小值、低射线鞍点的模式语义和可复核签名。",
        "已完成；验证 0 fail(s)",
    ),
    (
        "A3Rnew-P1 Jacobian 定版",
        "实现多尺度中心差分稳定平台判定；对非有限、病态、near-cutoff、扰动轨迹非 ground 等结果显式拒绝；修复细化尺度索引越界。",
        "已完成；正式路径为 central_reference_multiscale",
    ),
    (
        "A3Rnew-P2 模式完整性",
        "实现完整 beta/az 二维周期角域审计，O/X 独立落点图与检查点；联合拓扑绕数、目标邻近、映射折叠、Jacobian 奇异性和分支状态触发三级自适应加密，并以单元有界 TRF 求根。",
        "P2-0 至 P2-5 已完成",
    ),
    (
        "模式追踪与可恢复计算",
        "最终去重联合分支、状态、转向点、发射角、三维轨迹、反射高度和群时延；检查点记录配置、Provider 和代码哈希，不一致时拒绝恢复。",
        "已完成；专项验证 0 fail(s)",
    ),
]


RISKS = [
    "P2-6 冻结长算例尚未验收：需显式启用 adaptive_2d，完成 P0C、平滑磁离子 baseline、固定 Provider 哈希斜向 TID 的 O/X 全域审计，并确认连续两级模式集合稳定。",
    "完整 Jacobian 长回归运行超过 15 分钟无输出，已终止本轮进程；当前属于未完成验证，不能据此宣称通过或失败。",
    "变分 Jacobian 相对误差约 1.778e-2，仅保留为实验诊断；不宣称其已验证加速。",
    "P2 是同一 AH3D 求解器内部的模式发现完整性验证，不替代后续 A4 PHaRLAP 外部物理验证。",
    "near-cutoff 隧穿、模式转换、完整三维全波散射及绝对功率校准仍未解决。",
]


def set_fonts(doc: Document) -> None:
    normal = doc.styles["Normal"]
    normal.font.name = "SimSun"
    normal._element.rPr.rFonts.set(qn("w:eastAsia"), "SimSun")
    normal.font.size = Pt(10.5)
    for name in ("Title", "Heading 1", "Heading 2"):
        style = doc.styles[name]
        style.font.name = "SimHei"
        style._element.rPr.rFonts.set(qn("w:eastAsia"), "SimHei")


def add_table(doc: Document, headers: list[str], rows: list[tuple[str, ...]]) -> None:
    table = doc.add_table(rows=1, cols=len(headers))
    table.style = "Table Grid"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    for index, header in enumerate(headers):
        table.rows[0].cells[index].text = header
    for row in rows:
        cells = table.add_row().cells
        for index, value in enumerate(row):
            cells[index].text = value
            cells[index].vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.TOP
    for row in table.rows:
        for cell in row.cells:
            for paragraph in cell.paragraphs:
                for run in paragraph.runs:
                    run.font.size = Pt(9)


def build() -> None:
    doc = Document()
    set_fonts(doc)

    title = doc.add_heading("3D短波反射/散射混合传播模型周工作总结", 0)
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    subtitle = doc.add_paragraph("总结日期：2026年7月17日\n范围：上周工作精简回顾与本周 A3Rnew 实施进展")
    subtitle.alignment = WD_ALIGN_PARAGRAPH.CENTER

    doc.add_heading("一、上周工作（精简）", level=1)
    add_table(doc, ["工作项", "精简进展"], LAST_WEEK)

    doc.add_heading("二、本周工作", level=1)
    add_table(doc, ["工作项", "实施内容", "状态"], THIS_WEEK)

    doc.add_heading("三、验证结果", level=1)
    doc.add_paragraph(
        "P2 二维审计专项、P0C 高低射线、多尺度中心差分参考和 P2P 变分兼容测试均为 0 fail(s)；相关代码编译及 git diff --check 通过。"
    )

    doc.add_heading("四、当前进度与下一步", level=1)
    doc.add_paragraph(
        "当前执行位置为 A3Rnew-P2：P2-0 至 P2-5 已实现，下一步为 P2-6 冻结长算例验收。验收通过后进入 P3 文献/解析复现、P5new 精度验收，再进入 A4 PHaRLAP 外部验证。"
    )

    doc.add_heading("五、风险与未完成项", level=1)
    for item in RISKS:
        doc.add_paragraph(item, style="List Number")

    doc.save(OUTPUT)
    check = Document(OUTPUT)
    headings = [p.text for p in check.paragraphs if p.style.name.startswith("Heading")]
    assert OUTPUT.exists() and OUTPUT.stat().st_size > 0
    assert not any("参考文献" in text for text in headings)
    print(OUTPUT)


if __name__ == "__main__":
    build()
