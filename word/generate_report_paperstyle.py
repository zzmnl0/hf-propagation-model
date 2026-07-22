from pathlib import Path
import csv
from datetime import datetime
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT, WD_CELL_VERTICAL_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

ROOT = Path(r'D:\code11\IRI01\IRI03\TID')
CODE = ROOT / 'code'
OUT = CODE / 'output'
WORD = CODE / 'word'
WORD.mkdir(exist_ok=True)
DOCX = WORD / 'hf_hybrid_propagation_technical_report_paperstyle.docx'

SCENARIO_FILES = {
    '基准': 'modes_baseline.csv',
    'TID': 'modes_tid.csv',
    'Es': 'modes_es.csv',
    '等离子体泡': 'modes_bubble.csv',
    '综合扰动': 'modes_full.csv',
    'O/X': 'modes_OX.csv',
    '雷达': 'modes_radar.csv',
    'Spread-F': 'modes_spreadf.csv',
    '三维流管雷达': 'modes_tube_radar.csv',
}

def read_csv(name):
    p = OUT / name
    if not p.exists():
        return []
    with p.open('r', encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def val(row, key):
    try:
        x = row.get(key, '') if row else ''
        return None if x == '' else float(x)
    except Exception:
        return None

def fmt(x, n=3):
    return '-' if x is None else f'{x:.{n}f}'

sc = {k: read_csv(v) for k, v in SCENARIO_FILES.items()}
base = sc['基准'][0] if sc['基准'] else {}
tid = sc['TID'][0] if sc['TID'] else {}
base_tau, base_gp, base_pr, base_h = val(base,'tau_ms'), val(base,'group_path_km'), val(base,'Pr_dBW'), val(base,'h_reflect_km')
tid_tau, tid_gp, tid_pr, tid_h = val(tid,'tau_ms'), val(tid,'group_path_km'), val(tid,'Pr_dBW'), val(tid,'h_reflect_km')
ox = sc['O/X']
o_row = next((r for r in ox if r.get('label','').endswith('_O')), ox[0] if ox else {})
x_row = next((r for r in ox if r.get('label','').endswith('_X')), ox[1] if len(ox)>1 else {})
o_tau, x_tau, o_gp, x_gp = val(o_row,'tau_ms'), val(x_row,'tau_ms'), val(o_row,'group_path_km'), val(x_row,'group_path_km')
tube = sc['三维流管雷达']
tube_tau = [val(r,'tau_ms') for r in tube if val(r,'tau_ms') is not None]
tube_pr = [val(r,'Pr_dBW') for r in tube if val(r,'Pr_dBW') is not None]
tube_beta = [val(r,'beta_deg') for r in tube if val(r,'beta_deg') is not None]

def set_run(run, east='宋体', west='Times New Roman', size=12, bold=False, italic=False):
    run.font.name = west
    run._element.rPr.rFonts.set(qn('w:eastAsia'), east)
    run._element.rPr.rFonts.set(qn('w:ascii'), west)
    run._element.rPr.rFonts.set(qn('w:hAnsi'), west)
    run.font.size = Pt(size)
    run.bold = bold
    run.italic = italic

def setup_styles(doc):
    normal = doc.styles['Normal']
    normal.font.name = 'Times New Roman'
    normal._element.rPr.rFonts.set(qn('w:eastAsia'), '宋体')
    normal.font.size = Pt(12)
    for name, size in [('Heading 1', 14), ('Heading 2', 12), ('Heading 3', 12)]:
        st = doc.styles[name]
        st.font.name = 'Times New Roman'
        st._element.rPr.rFonts.set(qn('w:eastAsia'), '黑体')
        st.font.size = Pt(size)
        st.font.bold = True
        st.font.color.rgb = RGBColor(0,0,0)

def ptext(doc, text, first=True):
    p = doc.add_paragraph()
    r = p.add_run(text)
    set_run(r, size=12)
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    p.paragraph_format.line_spacing = 1.35
    p.paragraph_format.space_after = Pt(3)
    if first:
        p.paragraph_format.first_line_indent = Pt(24)
    return p

def heading(doc, text, level=1):
    p = doc.add_heading(text, level=level)
    for r in p.runs:
        set_run(r, east='黑体', size=14 if level==1 else 12, bold=True)
        r.font.color.rgb = RGBColor(0,0,0)
    p.paragraph_format.space_before = Pt(10 if level==1 else 6)
    p.paragraph_format.space_after = Pt(4)
    return p

def equation(doc, eq, num):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run(eq + f'    ({num})')
    set_run(r, east='Cambria Math', west='Cambria Math', size=11)
    p.paragraph_format.space_before = Pt(4)
    p.paragraph_format.space_after = Pt(4)
    return p

def caption(doc, text):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run(text)
    set_run(r, size=10.5)
    p.paragraph_format.space_after = Pt(6)
    return p

def shade(cell, color='D9EAF7'):
    tc = cell._tc.get_or_add_tcPr()
    shd = OxmlElement('w:shd')
    shd.set(qn('w:fill'), color)
    tc.append(shd)

def cell(cell, text, bold=False):
    cell.text = ''
    p = cell.paragraphs[0]
    r = p.add_run(str(text))
    set_run(r, size=10.5, bold=bold)
    cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER

def table(doc, title, headers, rows):
    caption(doc, title)
    t = doc.add_table(rows=1, cols=len(headers))
    t.style = 'Table Grid'
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i,h in enumerate(headers):
        cell(t.rows[0].cells[i], h, True); shade(t.rows[0].cells[i])
    for row in rows:
        cs = t.add_row().cells
        for i,v in enumerate(row): cell(cs[i], v)
    doc.add_paragraph()
    return t

def image(doc, name, cap, width=5.9):
    p = OUT / name
    if not p.exists():
        ptext(doc, f'图件 {name} 缺失，未插入。')
        return
    q = doc.add_paragraph(); q.alignment = WD_ALIGN_PARAGRAPH.CENTER
    q.add_run().add_picture(str(p), width=Inches(width))
    caption(doc, cap)

def result_row(name, rows, note):
    if not rows:
        return [name, '0', '-', '-', '-', note]
    r = rows[0]
    return [name, str(len(rows)), fmt(val(r,'tau_ms')), fmt(val(r,'group_path_km'),1), fmt(val(r,'Pr_dBW'),2), note]

# planning/audit update
plan_append = f'''
\n五、论文式修订记录（{datetime.now():%Y-%m-%d %H:%M:%S}）\n1. 按 scientific-writing 的 IMRAD 思路重构为“引言-模型与方法-数值实验-结果-讨论-结论”。\n2. 按 research-paper-writing 的贡献-证据-意义逻辑强化研究问题、方法创新和结果解释。\n3. 按 docx 规范重设 A4 页面、中文宋体小四正文、黑体标题、Cambria Math 公式、图题/表题小五号。\n4. 按 scientific-visualization 原则保留与结果讨论直接相关的输出图件。\n5. 按 paper-review 做自查：结果均来自 output/modes_*.csv 和已运行验证脚本；参考文献仍标注待核验。\n'''
plan_path = WORD / 'writing-plan.txt'
old_plan = plan_path.read_text(encoding='utf-8') if plan_path.exists() else ''
if '五、论文式修订记录' not in old_plan:
    plan_path.write_text(old_plan + plan_append, encoding='utf-8')

# document
doc = Document()
setup_styles(doc)
sec = doc.sections[0]
sec.page_width = Inches(8.27); sec.page_height = Inches(11.69)
sec.top_margin = Inches(0.9); sec.bottom_margin = Inches(0.85); sec.left_margin = Inches(1.0); sec.right_margin = Inches(1.0)

# title
t = doc.add_paragraph(); t.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = t.add_run('复杂电离层环境下 HF 短波反射/散射混合传播机制与三维建模')
set_run(r, east='黑体', size=16, bold=True)
st = doc.add_paragraph(); st.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = st.add_run('项目技术论文式报告')
set_run(r, size=12)
for line in [
    '项目：复杂电离层环境短波反射/散射混合传播机理及建模研究',
    '测试链路：TX 30N, 120E；RX 距离 1169 km；载频 10 MHz；时间 2020-06-01 12:00 UTC',
    f'生成时间：{datetime.now():%Y-%m-%d %H:%M:%S}'
]:
    pp = doc.add_paragraph(); pp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    rr = pp.add_run(line); set_run(rr, size=10.5)

doc.add_page_break()

heading(doc, '摘要', 1)
ptext(doc, f'复杂电离层环境中的 HF 短波传播同时受到背景电子密度、行进式电离层扰动、偶发 E 层、等离子体泡、Spread-F 不规则体、地磁场和 D 层吸收的共同影响。传统二维射线模型能够刻画主传播平面内的折返路径，但难以处理横向偏移、三维扰动和散射诊断；单一散射模型又难以保留通信和雷达链路所需的可达模式、群延迟和功率信息。为解决这一问题，本文构建并整理了一套三维 HF 短波反射/散射混合传播模型。该模型以 IRI-2016 为背景电子密度场，叠加 TID、Es、等离子体泡和 Spread-F 扰动，构造 RefractiveIndex3D 折射率场，并采用 Haselgrove RK4 与三维点对点 Newton shooting 求解传播模式。')
ptext(doc, f'在方法上，模型将传播输出统一组织为 SS1 至 SS13，包括群延迟、群路径、反射高度、通信接收功率、雷达回波功率、功率-时延谱、O/X 分裂、D 层吸收、散射诊断和横向偏移等参数。本文重点修正并强调 group_path_km = ∫ds/n 的群路径定义，避免将几何路径误用于时延和功率链路；同时采用 Appleton-Hartree 后处理计算 O/X 模式群路径差，并通过 scatter_coupler 为三维模式附加 Born/Rytov/MPS 散射诊断。数值验证表明，群路径、O/X 分裂、吸收、P2P Jacobian、散射耦合以及 Phase 5-7 的三维功能均通过验证。典型结果显示，基准 F 层单跳模式的群延迟为 {fmt(base_tau)} ms，TID 场景引起约 {fmt((tid_tau-base_tau)*1000 if tid_tau and base_tau else None,2)} μs 的时延变化，O/X 场景中 X 模相对 O 模晚约 {fmt((x_tau-o_tau)*1000 if x_tau and o_tau else None,2)} μs，三维流管雷达场景给出 {fmt(2*min(tube_tau) if tube_tau else None)} 至 {fmt(2*max(tube_tau) if tube_tau else None)} ms 的双程时延。结果说明，该框架能够同时支撑 HF 通信链路、OTH 雷达和定位误差分析。')
heading(doc, '关键词', 1)
ptext(doc, 'HF 短波传播；电离层；三维射线追踪；TID；O/X 分裂；功率-时延谱；Born/Rytov/MPS；OTH 雷达')

heading(doc, '1 引言', 1)
ptext(doc, 'HF 短波依赖电离层折射和反射实现越距传播，是远距离通信、超视距雷达和无源定位系统的重要物理基础。实际电离层具有明显的时空非均匀性：中尺度 TID 会改变主传播路径和到达角，偶发 E 层可能引入低高度反射或散射，F 区等离子体泡和 Spread-F 不规则体会引起相位起伏和功率谱展宽，D 层碰撞吸收则直接降低链路余量。这些效应通常同时存在，使单一机制模型难以解释观测到的时延、功率和角度误差。')
ptext(doc, '已有工程实现中，二维射线追踪常用于获得主路径和群延迟，散射模型常用于评估不规则体影响，雷达方程用于估计回波功率。然而，如果这些模块相互割裂，模型输出之间容易出现物理量定义不一致的问题。例如，若将几何路径直接用于群延迟或路径损耗，就会忽略冷等离子体中 n<1 对群传播路径的影响；若 O/X 模式只作为标签处理，则无法估计磁离子分裂造成的微小时延差；若三维横向扰动被投影到二维平面，则方位角误差无法直接输出。')
ptext(doc, '本文的研究目标是建立一套面向复杂电离层环境的三维混合传播建模流程，并将其整理为可复现的技术论文式报告。报告的贡献体现在三个方面：第一，给出从 IRI 背景到复杂扰动电子密度场、三维射线追踪和 P2P 求解的完整数据流；第二，统一定义 SS1 至 SS13 传播参数，使通信、雷达和定位分析共享同一物理输出；第三，基于现有验证脚本和输出文件解释典型场景结果，说明模型如何用于时延、功率、O/X 分裂和散射效应分析。')

heading(doc, '2 模型与方法', 1)
heading(doc, '2.1 电子密度场构建', 2)
ptext(doc, '模型以 IRI-2016 生成的电子密度廓线作为背景场。标准链路采用 2020-06-01 12:00 UTC、路径中点约 32.5N 和 120E 的 IRI 配置，垂直范围覆盖 60 至 600 km。背景电子密度被映射到 x-z 网格后，进一步叠加 TID、Es、等离子体泡和 Spread-F 扰动，从而形成用于射线追踪的 Ne_2d(x,z)。')
equation(doc, 'Ne(x,y,z,t) = Ne0(x,z) · [1 + A sin(kx x + ky y − ωt + φ)]', '1')
ptext(doc, '式(1)给出 TID 扰动的基本形式，其中 A 为相对扰动幅度，kx 和 ky 分别为沿链路方向和横向方向的波数。当前实现支持多分量 TID_3D 叠加。为避免 D 层近零密度区域压低 F 区扰动幅度，TID 相对扰动只在 Ne_bg 大于背景峰值 1% 的区域内归一化。对于多数场景，TID 先被烘焙进 Ne_2d，再由 RefractiveIndex3D 使用，从而降低射线步进时的重复计算成本。')
heading(doc, '2.2 三维折射率与射线追踪', 2)
ptext(doc, '冷等离子体各向同性折射率由电子等离子体频率决定。RefractiveIndex3D 在背景网格上预计算 n^2 及其空间梯度，并在 Haselgrove 方程右端通过插值查询。地球曲率采用局地展平修正，将直角坐标高度 z 转换为物理高度 z_phys。')
equation(doc, 'n²(x,z) = 1 − fp²(x,z) / f²', '2')
equation(doc, 'zphys = z − (x² + y²) / (2RE)', '3')
ptext(doc, '三维射线状态包含位置和波矢分量，使用 RK4 方法沿弧长推进。y 方向方程使模型能够直接计算横向落点 y_land 和方位角偏转，这是二维射线模型无法给出的量。点对点传播模式由 find_all_rays_3d 求解：算法先从射线扇面中选择接近接收点的候选射线，再以发射仰角 beta 和方位角 az 为未知量进行 Newton shooting。')
heading(doc, '2.3 群路径、功率链与后处理', 2)
ptext(doc, '本文将群路径定义为沿射线路径对 ds/n 的积分。该定义对应冷等离子体中 n_g≈1/n 的近似，是群延迟、通信路径损耗和雷达功率链的基础。几何路径 geom_path_km 仅表示轨迹长度，不直接用于群延迟。')
equation(doc, "P' = ∫ ds / n", '4')
equation(doc, "τ = P' / c", '5')
equation(doc, "Lfs = 20 log10(4πP' / λ)", '6')
ptext(doc, 'O/X 分裂通过 Appleton-Hartree 公式沿三维各向同性轨迹后处理。当前约定为 O 模使用 +disc 分支，X 模使用 -disc 分支。该方法可估计同一几何轨迹上的 O/X 群路径差，但还不是完整的三维 Appleton-Hartree 分支射线追踪。散射后处理由 scatter_coupler 完成，依据相位方差阈值将模式诊断为 Born、Rytov 或 MPS_required；默认情况下该诊断不直接改写主功率 Pr_W。')
equation(doc, 'n²O,X = 1 − X / [1 − ½YT²/(1−X) ± disc]', '7')

table(doc, '表1 模型模块与研究方法对应关系', ['模块', '研究方法作用', '主要输出'], [
    ['ionosphere_model.py', '构建 IRI 背景和复杂扰动电子密度场', 'Ne_2d'],
    ['ray_tracer_3d.py', '三维折射率查询与 Haselgrove RK4 追踪', '三维射线轨迹、群路径'],
    ['point_to_point_3d.py', 'P2P Newton shooting 求解可达模式', 'tau_ms、h_reflect_km、label'],
    ['scatter_coupler.py', 'Born/Rytov/MPS 散射诊断耦合', 'scatter_method、sigma_phi_sq、S4'],
    ['hybrid_model.py', '模式汇总、功率链和 P-D 谱生成', 'Pr_dBW、P-D spectrum、主模'],
])

table(doc, '表2 SS1 至 SS13 参数体系', ['类别', '参数', '物理意义'], [
    ['几何与时延', 'SS1 tau_ms；SS2 group_path_km；SS3 h_reflect_km/label', '描述主传播模式的时延、群路径、反射高度和跳数。'],
    ['通信链路', 'SS4 Pr_comm；SS9 F_focus', '评估通信接收功率和射线束聚焦/散焦影响。'],
    ['雷达链路', 'SS5 Pr_radar；SS7 P-D spectrum', '给出 OTH 雷达双程回波功率和功率-时延谱结构。'],
    ['磁离子效应', 'SS8 delta_tau_OX', '估计 O/X 模式群路径差和定位误差来源。'],
    ['复杂传播', 'SS6、SS10、SS11、SS12、SS13', '描述 Es、Born/Rytov/MPS 散射、D 层吸收和横向偏移。'],
])

heading(doc, '3 数值实验设计与验证', 1)
ptext(doc, '数值实验采用标准链路 TX 30N、120E 至 1169 km 接收点，载频 10 MHz，时间为 2020-06-01 12:00 UTC。实验分为两类：第一类为验证脚本，用于检查物理公式、数值一致性和模块接口；第二类为场景输出，用于比较基准、TID、Es、bubble、O/X、雷达和 Spread-F 条件下的模式参数与图件。')
table(doc, '表3 验证脚本与验证结论', ['验证脚本', '验证内容', '结论'], [
    ['_verify_group_path.py', 'group_path=∫ds/n，τ=P\'/c，3D/2D 一致性', '0 fail(s)'],
    ['_verify_OX.py', 'AH 各向同性极限、O/X 模式字段和时延差', 'PASSED'],
    ['_verify_absorption_3d.py', '通信单程吸收与雷达双程吸收', '0 fail(s)'],
    ['_verify_p2p_jacobian.py', 'forward 与 central Jacobian 一致性', '0 fail(s)'],
    ['_verify_scatter_coupling.py', 'Born/Rytov/MPS 阈值选择', '0 fail(s)'],
    ['_verify_phase5.py', '三维背景、IGRF、多分量 TID', '单独重跑 0 fail(s)'],
    ['_verify_phase6.py', '三维射线、P2P、Born/Rytov', '0 fail(s)'],
    ['_verify_phase7.py', '3D bubble、tube、Booker-Gordon、MPS', '0 fail(s)'],
])
ptext(doc, 'Phase 5 在与其他脚本并行运行时曾触发 ppigrf 内存分配异常，单独重跑后通过。该现象反映的是运行资源竞争或第三方库内存分配问题，而非模型物理公式失败。因此，本文在结果讨论中采用单独重跑后的验证结论，同时保留该异常记录以保证可复现性。')

heading(doc, '4 结果', 1)
heading(doc, '4.1 基准链路传播结果', 2)
ptext(doc, f'基准场景输出一条主要 F 层低角传播模式。该模式 tau_ms 为 {fmt(base_tau)} ms，group_path_km 为 {fmt(base_gp,1)} km，反射高度为 {fmt(base_h,1)} km，接收功率为 {fmt(base_pr,2)} dBW。由于地面水平距离为 1169 km，而群路径达到 {fmt(base_gp,1)} km，结果说明折返轨迹长度和 n<1 的群折射率修正共同增加了等效传播路径。')
image(doc, 'ray_fan_baseline.png', '图1 基准场景三维射线扇面。射线路径显示主要可达 F 层低角模式。')
heading(doc, '4.2 TID 对时延与横向偏移的影响', 2)
ptext(doc, f'TID 场景中，tau_ms 从基准值 {fmt(base_tau)} ms 增加到 {fmt(tid_tau)} ms，变化约 {fmt((tid_tau-base_tau)*1000 if tid_tau and base_tau else None,2)} μs；group_path_km 从 {fmt(base_gp,1)} km 增加到 {fmt(tid_gp,1)} km，变化约 {fmt(tid_gp-base_gp if tid_gp and base_gp else None,2)} km。该量级说明当前 TID 配置对主模时延的影响较小，但它能够改变折返高度和横向相位梯度。对于定位问题，即使微秒级时延变化和数公里级横向偏移也可能转化为显著的位置误差。')
image(doc, 'ray_deflect_tid.png', '图2 TID 引起的横向偏移诊断。该图用于解释 SS13 的 y_land 和方位角偏转。')
heading(doc, '4.3 O/X 模式分裂结果', 2)
ptext(doc, f'O/X 场景中，O 模 tau_ms 为 {fmt(o_tau)} ms，X 模 tau_ms 为 {fmt(x_tau)} ms，X 模相对 O 模晚约 {fmt((x_tau-o_tau)*1000 if x_tau and o_tau else None,2)} μs；O/X 群路径差约 {fmt(x_gp-o_gp if x_gp and o_gp else None,2)} km。该结果表明，磁离子效应在 10 MHz 链路上虽然不会改变主传播模式类型，但会产生可解析的群延迟差。对于高精度 TDOA 或雷达测距，该差异不应被简单忽略。')
image(doc, 'pd_OX.png', '图3 O/X 场景功率-时延谱。谱峰间隔反映磁离子模式群路径差。')
heading(doc, '4.4 雷达与流管功率-时延结果', 2)
ptext(doc, f'三维流管雷达场景给出多条相邻仰角模式。发射仰角约为 {fmt(min(tube_beta),1) if tube_beta else "-"} 至 {fmt(max(tube_beta),1) if tube_beta else "-"} 度，单程 tau_ms 范围为 {fmt(min(tube_tau)) if tube_tau else "-"} 至 {fmt(max(tube_tau)) if tube_tau else "-"} ms，对应双程时延为 {fmt(2*min(tube_tau)) if tube_tau else "-"} 至 {fmt(2*max(tube_tau)) if tube_tau else "-"} ms；回波功率约为 {fmt(min(tube_pr),2) if tube_pr else "-"} 至 {fmt(max(tube_pr),2) if tube_pr else "-"} dBW。结果说明，雷达回波谱不是单一主模点，而是由相邻射线管、双程传播损耗和聚焦/散焦共同塑造的分布。')
image(doc, 'pd_tube_radar.png', '图4 三维流管雷达功率-时延谱。多条相邻模式形成可解释的回波时延结构。')
heading(doc, '4.5 复杂扰动场景对比', 2)
rows=[]
notes={'基准':'F 层主模基准。','TID':'时延和路径微扰。','Es':'低高度薄层改变模式和谱形态。','等离子体泡':'F 区耗空导致局部折射率变化。','综合扰动':'多机制叠加 pipeline 稳定性。','O/X':'磁离子分裂。','雷达':'双程雷达功率链。','Spread-F':'不规则体与相位屏诊断。','三维流管雷达':'相邻仰角射线管和回波谱。'}
for name, rows_csv in sc.items():
    rows.append(result_row(name, rows_csv, notes.get(name,'')))
table(doc, '表4 典型场景输出结果摘要', ['场景', '模式数', '代表 τ/ms', '代表 P\'/km', '代表 Pr/dBW', '结果解读'], rows)
image(doc, 'ray_fan_spreadf.png', '图5 Spread-F 场景射线扇面。图中轨迹用于辅助理解不规则体条件下的传播模式变化。')

heading(doc, '5 讨论', 1)
ptext(doc, '结果首先验证了群路径定义在整个传播链中的核心作用。基准链路中，group_path_km 明显大于地面水平距离，说明 HF 链路的时延和功率损耗不能由几何距离直接替代。该修正对于 P-D 谱横轴、雷达测距和通信链路预算均具有基础意义。')
ptext(doc, 'TID 结果显示，在当前扰动幅度和波长配置下，主模时延变化处于微秒量级，但三维模型能够进一步给出横向偏移和方位角偏转。这一点是三维建模相对于二维模型的关键增量：二维结果可解释沿链路方向的折返路径，而三维结果能够服务于定位误差分析。')
ptext(doc, 'O/X 结果说明，磁离子分裂对群延迟的影响虽小，但在高精度系统中不可忽略。当前实现沿各向同性三维轨迹重新积分 Appleton-Hartree 折射率，适合进行快速诊断和误差估计；但严格的 O/X 分支轨迹分离仍需要完整三维 Appleton-Hartree 射线追踪。')
ptext(doc, '散射诊断和流管雷达结果表明，功率-时延谱不仅由主射线路径决定，还受不规则体相位方差、流管截面积、聚焦/散焦和双程传播损耗共同影响。当前 scatter_coupler 只附加 Born/Rytov/MPS 诊断而不改写 Pr_W，这一设计使主传播链保持稳定，但在强散射条件下仍需结合 PE/SSF 或多相位屏结果进行功率修正。')
ptext(doc, '本报告的主要局限在于参考文献尚未核验，且结果分析主要基于现有场景输出而非系统参数扫描。后续若面向正式投稿，应补充多频率、多时刻、多扰动强度的敏感性实验，给出统计误差或不确定性范围，并与实测 P-D 谱或已有传播模型进行定量比较。')

heading(doc, '6 结论', 1)
ptext(doc, '本文将复杂电离层 HF 短波混合传播代码整理为论文式技术报告，系统阐述了电子密度场构建、三维射线追踪、点对点求解、功率链、O/X 后处理、吸收和散射诊断之间的数据流关系。通过 SS1 至 SS13 参数体系，模型能够同时服务通信链路余量、OTH 雷达功率-时延谱和定位误差分析。')
ptext(doc, '基于验证脚本和典型输出，当前模型已在群路径、O/X 分裂、三维吸收、P2P Jacobian、散射耦合和三维扩展方面通过检查。结果表明，该框架能够解释基准 F 层传播、TID 微扰、磁离子分裂、流管雷达回波谱和 Spread-F 条件下的模式变化。后续工作应集中于完整三维 AH 射线追踪、强散射功率反馈、实测谱反演和参考文献核验。')

heading(doc, '参考文献（待核验）', 1)
for ref in ['IRI-2016 国际参考电离层模型相关文献。','Haselgrove 射线追踪与 Jones & Stephenson 1975 地球展平/射线追踪相关文献。','Appleton-Hartree 磁离子理论相关文献。','Coleman 1997/1998 flux-tube 聚焦与 OTH 雷达传播相关文献。','Hao 2017 偶发 E 层三段传播模型相关文献。','Rino 1979 电离层不规则体相位屏和闪烁理论相关文献。','Knepp 多相位屏模型相关文献。','Pederick & Cervera 2014 D 层吸收参数化相关文献。']:
    ptext(doc, ref)

heading(doc, '附录A 主要配置与输出文件', 1)
table(doc, '表A1 关键配置参数', ['参数组', '字段', '说明'], [
    ['基础链路', 'FREQ_MHZ, TX_POS, RX_POS', '载频、发射点和接收点。'],
    ['背景网格', 'BG_X, BG_Z', '电子密度背景网格。'],
    ['RT_3D', 'ds_km, beta_min, beta_max, n_fan', '三维射线步长和扇面。'],
    ['P2P_3D', 'n_init, tol_km, d_beta, d_az, jacobian', '三维 P2P Newton shooting。'],
    ['TID/TID_3D', 'amplitude, wavelength, period, azimuth', 'TID 扰动参数。'],
    ['SCATTER', 'Cs_rel, L_outer_km, p, DeltaL_km', '散射诊断参数。'],
    ['ABSORPTION', 'enable, A0, chi_deg', 'D 层吸收参数。'],
])
ptext(doc, '主要运行命令为 conda run -n pytorch_cpu python main.py。验证命令包括 tests/_verify_group_path.py、tests/_verify_OX.py、tests/_verify_absorption_3d.py、tests/_verify_p2p_jacobian.py、tests/_verify_scatter_coupling.py、tests/_verify_phase5.py、tests/_verify_phase6.py 和 tests/_verify_phase7.py。输出文件包括 ray_fan_*.png、pd_*.png、ray_deflect_tid.png 和 modes_*.csv。')

# footer
for s in doc.sections:
    ftr = s.footer.paragraphs[0]
    ftr.alignment = WD_ALIGN_PARAGRAPH.CENTER
    rr = ftr.add_run('复杂电离层 HF 短波混合传播技术报告')
    set_run(rr, size=9)

doc.save(DOCX)
chk = Document(DOCX)
validation = f'docx={DOCX}\nparagraphs={len(chk.paragraphs)}\ntables={len(chk.tables)}\nsize_bytes={DOCX.stat().st_size}\nrevision=paper-style-v2\n'
(WORD / 'validation.txt').write_text(validation, encoding='utf-8')
review = '''论文式修订自查\n\n1. 结构：已由说明书式章节调整为摘要、引言、模型与方法、数值实验与验证、结果、讨论、结论和附录。\n2. 语言：正文改为连续科研段落，减少项目符号和代码说明式表达。\n3. 方法：集中说明电子密度场、三维折射率、RK4 射线追踪、P2P 求解、群路径、功率链、O/X 后处理和散射诊断。\n4. 结果：增加基准、TID、O/X、流管雷达和复杂扰动对比的数值解读。\n5. 格式：正文宋体小四，标题黑体，公式居中并使用 Cambria Math，图题表题小五号。\n6. 限制：参考文献仍为待核验；O/X 明确说明不是完整三维 AH 射线追踪；scatter_coupler 默认不改写 Pr_W。\n'''
(WORD / 'paper_review_notes.txt').write_text(review, encoding='utf-8')
print('OK paper-style-v2')
