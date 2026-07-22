from pathlib import Path
import shutil
import pandas as pd
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
from docx.shared import Inches, Pt, RGBColor
from docx.oxml.ns import qn

ROOT = Path.cwd()
WORD = ROOT / 'code' / 'word'
OUT = ROOT / 'code' / 'output'
TARGET = WORD / 'hf_hybrid_propagation_technical_report_paperstyle.docx'
SOURCE = WORD / 'hf_hybrid_propagation_technical_report_results_enhanced.docx'
BACKUP = WORD / 'hf_hybrid_propagation_technical_report_paperstyle.pre_results_merge_backup.docx'

summary = pd.read_csv(OUT / 'report_scenario_summary.csv')
detail = pd.read_csv(OUT / 'report_mode_details.csv')

def row(scenario):
    rows = summary[summary['scenario'].astype(str).str.lower() == scenario.lower()]
    return rows.iloc[0] if len(rows) else pd.Series(dtype=object)

def fmt(value, digits=3, missing='-'):
    try:
        if pd.isna(value):
            return missing
        return f'{float(value):.{digits}f}'
    except Exception:
        return str(value) if value not in (None, '') else missing

def set_run(run, size=12, bold=False):
    run.font.name = 'Times New Roman'
    run._element.rPr.rFonts.set(qn('w:eastAsia'), '宋体')
    run._element.rPr.rFonts.set(qn('w:ascii'), 'Times New Roman')
    run._element.rPr.rFonts.set(qn('w:hAnsi'), 'Times New Roman')
    run.font.size = Pt(size)
    run.bold = bold
    run.font.color.rgb = RGBColor(0, 0, 0)

def style_paragraph(p, first_line=True, center=False):
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER if center else WD_ALIGN_PARAGRAPH.JUSTIFY
    p.paragraph_format.line_spacing = 1.35
    p.paragraph_format.space_after = Pt(3)
    if first_line and not center:
        p.paragraph_format.first_line_indent = Pt(24)
    for r in p.runs:
        set_run(r, 12)
    return p

def remove_chapter_4(doc):
    paras = doc.paragraphs
    start = next(i for i, p in enumerate(paras) if p.text.strip().startswith('4 结果'))
    end = next(i for i, p in enumerate(paras[start + 1:], start + 1) if p.text.strip().startswith('5 '))
    body = doc._body._element
    start_el = paras[start]._element
    end_el = paras[end]._element
    children = list(body)
    start_idx = children.index(start_el)
    end_idx = children.index(end_el)
    for el in children[start_idx:end_idx]:
        body.remove(el)
    return end_el

def insert_element_before(doc, marker_el, element):
    body = doc._body._element
    body.remove(element)
    body.insert(list(body).index(marker_el), element)

def add_heading_before(doc, marker_el, text, level):
    p = doc.add_heading(text, level=level)
    p.paragraph_format.space_before = Pt(10 if level == 1 else 6)
    p.paragraph_format.space_after = Pt(4)
    for r in p.runs:
        set_run(r, 14 if level == 1 else 12, bold=True)
    insert_element_before(doc, marker_el, p._element)
    return p

def add_para_before(doc, marker_el, text):
    p = doc.add_paragraph()
    p.add_run(text)
    style_paragraph(p)
    insert_element_before(doc, marker_el, p._element)
    return p

def add_caption_before(doc, marker_el, text):
    p = doc.add_paragraph()
    p.add_run(text)
    style_paragraph(p, first_line=False, center=True)
    for r in p.runs:
        set_run(r, 10.5)
    insert_element_before(doc, marker_el, p._element)
    return p

def add_image_before(doc, marker_el, filename, caption, width=6.2):
    path = OUT / filename
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    if path.exists():
        p.add_run().add_picture(str(path), width=Inches(width))
    else:
        p.add_run(f'[图件缺失：{filename}]')
    insert_element_before(doc, marker_el, p._element)
    add_caption_before(doc, marker_el, caption)

def add_table_before(doc, marker_el, title, headers, rows):
    add_caption_before(doc, marker_el, title)
    table = doc.add_table(rows=1, cols=len(headers))
    table.style = 'Table Grid'
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, h in enumerate(headers):
        cell = table.rows[0].cells[i]
        cell.text = str(h)
        cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
        for p in cell.paragraphs:
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            for r in p.runs:
                set_run(r, 9.5, bold=True)
    for row_values in rows:
        cells = table.add_row().cells
        for i, value in enumerate(row_values):
            cells[i].text = str(value)
            cells[i].vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            for p in cells[i].paragraphs:
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                for r in p.runs:
                    set_run(r, 9)
    insert_element_before(doc, marker_el, table._element)
    spacer = doc.add_paragraph()
    insert_element_before(doc, marker_el, spacer._element)

if not TARGET.exists():
    raise FileNotFoundError(TARGET)
if not SOURCE.exists():
    raise FileNotFoundError(SOURCE)
shutil.copy2(TARGET, BACKUP)
doc = Document(TARGET)
marker = remove_chapter_4(doc)
base = row('baseline')
tid = row('tid')
es = row('es')
bubble = row('bubble')
full = row('full')
ox = row('OX')
radar = row('radar')
spreadf = row('spreadf')
tube = detail[detail['tau_2way_ms'].notna()].copy()

add_heading_before(doc, marker, '4 结果', 1)
add_heading_before(doc, marker, '4.1 场景主模对比', 2)
rows = []
for _, r in summary.iterrows():
    rows.append([
        r.get('scenario_label', r.get('scenario', '')),
        int(r.get('n_modes', 0)) if not pd.isna(r.get('n_modes', None)) else '-',
        r.get('label', '-'),
        r.get('wave_mode', '-'),
        fmt(r.get('tau_ms'), 3),
        fmt(r.get('group_path_km'), 1),
        fmt(r.get('h_reflect_km'), 1),
        fmt(r.get('Pr_dBW'), 2),
    ])
add_table_before(doc, marker, '表4-1 各场景主模结果汇总', ['场景', '模式数', '主模', '波模', 'tau/ms', "P'/km", 'href/km', 'Pr/dBW'], rows)
add_para_before(doc, marker, f"基准场景的主模群延迟为 {fmt(base.get('tau_ms'), 3)} ms，群路径为 {fmt(base.get('group_path_km'), 1)} km，反射高度为 {fmt(base.get('h_reflect_km'), 1)} km，接收功率为 {fmt(base.get('Pr_dBW'), 2)} dBW。这里的 group_path_km 表示冷等离子体近似下的积分量 integral(ds/n)，因此它不是几何路程，而是同时包含折返轨迹长度和折射率修正的传播量。")
add_para_before(doc, marker, f"从主模表可以看出，Es 与综合扰动场景优先选中低高度 Es 模式，群延迟分别相对基准变化 {fmt(es.get('delta_tau_us'), 2)} us 和 {fmt(full.get('delta_tau_us'), 2)} us；TID、bubble 与 Spread-F 场景仍保持 1F 低角模式，但群路径和时延发生小幅偏移。这种对比说明，薄层型结构主要改变可达模态集合，而缓变 F 区扰动更多表现为同一模态参数的连续漂移。")
add_image_before(doc, marker, 'report_fig_scenario_comparison.png', '图4-1 各场景主模的群延迟、群路径、接收功率和反射高度对比。')

add_heading_before(doc, marker, '4.2 相对基准差异', 2)
add_para_before(doc, marker, f"TID 场景相对基准的时延差为 {fmt(tid.get('delta_tau_us'), 2)} us，群路径差为 {fmt(tid.get('delta_group_path_km'), 3)} km，功率差为 {fmt(tid.get('delta_Pr_dB'), 2)} dB。该量级表明多方向 TID 在本链路中主要引起微秒级群延迟扰动和公里量级以内的等效群路径变化，功率链路变化较小。")
add_para_before(doc, marker, f"Bubble 与 Spread-F 场景的主模仍接近 TID 类扰动的连续偏移特征，其中 bubble 的群路径差为 {fmt(bubble.get('delta_group_path_km'), 3)} km，Spread-F 的群路径差为 {fmt(spreadf.get('delta_group_path_km'), 3)} km。相反，Es 和综合扰动场景的相对差异更大，这是因为主导模态从 F 区折返切换到低高度 Es 相关路径，导致反射高度、群路径和时延同步下降。")
add_image_before(doc, marker, 'report_fig_delta_vs_baseline.png', '图4-2 各扰动场景相对基准的群延迟、群路径和接收功率差异。')

add_heading_before(doc, marker, '4.3 O/X 模式分裂结果', 2)
add_para_before(doc, marker, f"O/X 场景在同一三维各向同性轨迹上进行 Appleton-Hartree 后处理，得到 O 波与 X 波之间的群延迟差 {fmt(ox.get('delta_tau_OX_us'), 3)} us，群路径差 {fmt(ox.get('delta_group_path_OX_km'), 3)} km。该结果可作为 10 MHz 链路上磁离子分裂引起的 TDOA 量级估计，但其物理含义应限定为后处理近似，而不是完整的三维 AH 分支射线追踪。")
add_para_before(doc, marker, "图4-3 将 O/X 的群延迟、群路径及二者差值并列展示。由于 O/X 处理沿用同一空间轨迹，图中差异主要反映折射率分支积分项，而不包含两分支独立弯曲、横向偏转和焦散结构的反馈效应。后续若要讨论强磁场各向异性传播，应将 O/X 分支纳入三维射线方程本身。")
add_image_before(doc, marker, 'report_fig_ox_split.png', '图4-3 O/X 模式的群延迟、群路径及分裂差值。')

add_heading_before(doc, marker, '4.4 三维流管雷达仰角扫描', 2)
if len(tube):
    add_para_before(doc, marker, f"三维流管雷达结果覆盖仰角 {tube['beta_deg'].min():.1f}-{tube['beta_deg'].max():.1f} deg，对应双程时延 {tube['tau_2way_ms'].min():.3f}-{tube['tau_2way_ms'].max():.3f} ms，接收功率 {tube['Pr_dBW'].min():.2f}-{tube['Pr_dBW'].max():.2f} dBW。随仰角变化，群路径和双程时延呈系统性扫描趋势，功率变化则同时受到传播损耗和流管聚焦或散焦的影响。")
else:
    add_para_before(doc, marker, "三维流管雷达结果未在明细表中给出 tau_2way_ms，本文仅保留图件和输出文件索引，不对双程时延范围作额外结论。")
add_para_before(doc, marker, "该雷达结果说明回波谱不是单一主模点，而是由相邻射线管、双程传播损耗、流管截面积和三维聚焦因子共同塑造的分布。与通信链路不同，雷达场景需要同时报告单程群路径、双程时延和功率扫描曲线，才能支撑后续对回波展宽和角度敏感性的讨论。")
add_image_before(doc, marker, 'report_fig_tube_radar_scan.png', '图4-4 三维流管雷达的仰角-双程时延、仰角-功率和仰角-群路径关系。')

add_heading_before(doc, marker, '4.5 复杂扰动场景综合矩阵', 2)
add_para_before(doc, marker, "结果矩阵将时延差、群路径差、功率差和散射诊断量归一化后放在同一图中，用于横向比较各类扰动对传播链路的主导影响。矩阵中的空白或接近零值应结合原始 CSV 字段解释：空白通常表示该场景未产生相应诊断量，接近零则表示相对基准的变化在本仿真设置下较弱。")
add_para_before(doc, marker, "综合比较表明，TID、bubble 和 Spread-F 更适合用于检验同一 F 区主模的时延稳定性与横向偏转，Es 与 full 场景则强调模态竞争和低高度反射路径的引入；散射耦合字段提供 Born、Rytov 或 MPS 诊断，但默认不直接改写 Pr_W，因此功率结论仍主要来自射线路径损耗、吸收和雷达几何项。")
add_image_before(doc, marker, 'report_fig_result_matrix.png', '图4-5 场景-指标影响矩阵，用于综合比较不同扰动机制的结果特征。')

OUTPUT = WORD / 'hf_hybrid_propagation_technical_report_paperstyle_results_integrated.docx'
doc.save(OUTPUT)
chk = Document(OUTPUT)
imgs = [r for r in chk.part.rels.values() if 'image' in r.reltype]
print('merged_target=' + str(OUTPUT))
print('backup=' + str(BACKUP))
print('paragraphs=%d tables=%d images=%d size=%d' % (len(chk.paragraphs), len(chk.tables), len(imgs), OUTPUT.stat().st_size))

