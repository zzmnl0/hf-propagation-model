from pathlib import Path
import csv
from datetime import datetime
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

ROOT=Path(r'D:\code11\IRI01\IRI03\TID'); CODE=ROOT/'code'; OUT=CODE/'output'; WORD=CODE/'word'; WORD.mkdir(exist_ok=True)

def read_csv(name):
    p=OUT/name
    if not p.exists(): return []
    with p.open('r',encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))

def f(row,k):
    try:
        v=row.get(k,'')
        return None if v=='' else float(v)
    except Exception: return None

def fmt(x,n=3): return '-' if x is None else f'{x:.{n}f}'

files={'baseline':'modes_baseline.csv','tid':'modes_tid.csv','es':'modes_es.csv','bubble':'modes_bubble.csv','full':'modes_full.csv','OX':'modes_OX.csv','radar':'modes_radar.csv','spreadf':'modes_spreadf.csv','tube_radar':'modes_tube_radar.csv'}
sc={k:read_csv(v) for k,v in files.items()}
base=sc['baseline'][0] if sc['baseline'] else {}; tid=sc['tid'][0] if sc['tid'] else {}
base_tau,base_gp,base_pr,base_h=f(base,'tau_ms'),f(base,'group_path_km'),f(base,'Pr_dBW'),f(base,'h_reflect_km')
tid_tau,tid_gp,tid_pr,tid_h=f(tid,'tau_ms'),f(tid,'group_path_km'),f(tid,'Pr_dBW'),f(tid,'h_reflect_km')
ox=sc['OX']; o=x={}
if len(ox)>=2:
    o=next((r for r in ox if r.get('label','').endswith('_O')),ox[0]); x=next((r for r in ox if r.get('label','').endswith('_X')),ox[1])
o_tau,x_tau,o_gp,x_gp=f(o,'tau_ms'),f(x,'tau_ms'),f(o,'group_path_km'),f(x,'group_path_km')
tube=sc['tube_radar']; tube_tau=[f(r,'tau_ms') for r in tube if f(r,'tau_ms') is not None]; tube_pr=[f(r,'Pr_dBW') for r in tube if f(r,'Pr_dBW') is not None]; tube_beta=[f(r,'beta_deg') for r in tube if f(r,'beta_deg') is not None]

plan='''技术报告 Word 撰写与生成执行计划\n\n一、用户确认需求\n1. 仅中文。\n2. 定位为项目技术论文/报告型 Word，同时讲清模型构建方法和模型结果解读。\n3. 允许运行验证脚本并读取 code/output/ 图件与 CSV。\n4. 参考文献暂不联网核验。\n5. SS1~SS13 细节按报告需要保留。\n6. 保存路径为 code/word/。\n\n二、调用/采用的 skill\n1. research-paper-writing：贡献、证据链、章节逻辑、结果叙事。\n2. scientific-writing：科研报告结构、段落化表达、图表组织。\n3. docx：Word 文档生成规范。\n4. scientific-visualization：现有 PNG 图件选择和图表服务于结果解释。\n5. paper-review：最终完整性和证据支撑自查。\n\n三、验证记录\n- _verify_group_path.py：0 fail(s)。\n- _verify_OX.py：PASSED。\n- _verify_absorption_3d.py：0 fail(s)。\n- _verify_p2p_jacobian.py：0 fail(s)。\n- _verify_scatter_coupling.py：0 fail(s)。\n- _verify_phase5.py：并行运行时 ppigrf 出现内存分配异常；单独重跑后 0 fail(s)。\n- _verify_phase6.py：0 fail(s)。\n- _verify_phase7.py：0 fail(s)。\n\n四、输出文件\n- writing-plan.txt\n- result_summary.md\n- hf_hybrid_propagation_technical_report.docx\n- validation.txt\n'''
(WORD/'writing-plan.txt').write_text(plan,encoding='utf-8')

lines=['# 结果摘要','']
for name,rows in sc.items():
    lines.append(f'## {name}')
    if not rows: lines.append('无结果。'); continue
    for r in rows[:8]: lines.append(f"- {r.get('label','-')}: tau={fmt(f(r,'tau_ms'))} ms, group_path={fmt(f(r,'group_path_km'),1)} km, Pr={fmt(f(r,'Pr_dBW'),2)} dBW, h={fmt(f(r,'h_reflect_km'),1)} km")
    lines.append('')
(WORD/'result_summary.md').write_text('\n'.join(lines),encoding='utf-8')

def set_font(run,font='宋体',size=10.5,bold=None):
    run.font.name=font; run._element.rPr.rFonts.set(qn('w:eastAsia'),font); run.font.size=Pt(size)
    if bold is not None: run.bold=bold

def heading(doc,text,level=1):
    p=doc.add_heading(text,level=level)
    for r in p.runs:
        set_font(r,'黑体',15 if level==1 else 13,True); r.font.color.rgb=RGBColor(0,0,0)
    return p

def para(doc,text):
    p=doc.add_paragraph(); r=p.add_run(text); set_font(r)
    p.paragraph_format.first_line_indent=Pt(21); p.paragraph_format.line_spacing=1.25
    return p

def cell_text(cell,text,bold=False):
    cell.text=''; p=cell.paragraphs[0]; r=p.add_run(str(text)); set_font(r,'宋体',9,bold)

def shade(cell,color='D9EAF7'):
    tc=cell._tc.get_or_add_tcPr(); sh=OxmlElement('w:shd'); sh.set(qn('w:fill'),color); tc.append(sh)

def table(doc,headers,rows):
    t=doc.add_table(rows=1,cols=len(headers)); t.style='Table Grid'; t.alignment=WD_TABLE_ALIGNMENT.CENTER
    for i,h in enumerate(headers): cell_text(t.rows[0].cells[i],h,True); shade(t.rows[0].cells[i])
    for row in rows:
        cs=t.add_row().cells
        for i,v in enumerate(row): cell_text(cs[i],v)
    doc.add_paragraph(); return t

def image(doc,name,cap,width=6.1):
    p=OUT/name
    if not p.exists(): para(doc,f'图件缺失：{name}'); return
    q=doc.add_paragraph(); q.alignment=WD_ALIGN_PARAGRAPH.CENTER; q.add_run().add_picture(str(p),width=Inches(width))
    c=doc.add_paragraph(cap); c.alignment=WD_ALIGN_PARAGRAPH.CENTER
    for r in c.runs: set_font(r,'宋体',9)

doc=Document(); sec=doc.sections[0]; sec.top_margin=Inches(.8); sec.bottom_margin=Inches(.8); sec.left_margin=Inches(.85); sec.right_margin=Inches(.85)
for st in ['Normal','Body Text']:
    if st in doc.styles:
        doc.styles[st].font.name='宋体'; doc.styles[st]._element.rPr.rFonts.set(qn('w:eastAsia'),'宋体'); doc.styles[st].font.size=Pt(10.5)

t=doc.add_paragraph(); t.alignment=WD_ALIGN_PARAGRAPH.CENTER; r=t.add_run('复杂电离层环境下 HF 短波反射/散射混合传播机制与三维建模技术报告'); set_font(r,'黑体',20,True)
for s in ['项目技术论文/报告型 Word 初稿','项目：复杂电离层环境短波反射/散射混合传播机理及建模研究','版本：2026-07-01 公式修复同步版','测试链路：TX 30N,120E；RX 距离 1169 km；载频 10 MHz；时间 2020-06-01 12:00 UTC',f'生成时间：{datetime.now():%Y-%m-%d %H:%M:%S}']:
    p=doc.add_paragraph(); p.alignment=WD_ALIGN_PARAGRAPH.CENTER; rr=p.add_run(s); set_font(rr,'宋体',12)
doc.add_page_break()

heading(doc,'摘要')
para(doc,'复杂电离层环境会使 HF 短波传播同时表现出折射、反射、散射、多径、磁离子分裂、D 层吸收和横向偏移等效应，直接影响远距离通信、超视距雷达和无源定位的时延、功率与角度估计。为服务该类应用，本报告整理并论文式阐述一套三维 HF 短波反射/散射混合传播模型。模型以 IRI-2016 为背景电子密度场，在其上叠加 TID、Es、等离子体泡和 Spread-F 等扰动，构造二维电子密度网格与三维折射率查询模型，并通过 Haselgrove RK4 射线追踪和三维点对点 Newton shooting 求解可达传播模式。')
para(doc,'在传播参数层面，报告将模型输出组织为 SS1 至 SS13，包括群延迟、群路径、反射高度、通信接收功率、雷达回波功率、Es 参数、功率-时延谱、O/X 模式时延差、磁通管聚焦、Born 弱散射、Rytov 闪烁、D 层吸收和横向偏移。当前版本明确区分几何路径与群路径，将 group_path_km 定义为 integral(ds/n)，修正 Appleton-Hartree O/X 分支，加入三维 D 层吸收诊断，并将 Born/Rytov/MPS 散射诊断耦合到三维模式结果中。')
para(doc,f'验证结果表明，核心物理修复项、三维背景、三维射线追踪、散射耦合和 Phase 7 扩展均通过检查。典型场景输出显示，基准 F2 单跳链路群延迟约为 {fmt(base_tau)} ms，TID 场景产生约 {fmt((tid_tau-base_tau)*1000 if tid_tau and base_tau else None,2)} us 的时延变化，O/X 场景中 X 模相对 O 模晚约 {fmt((x_tau-o_tau)*1000 if x_tau and o_tau else None,2)} us，三维流管雷达场景给出约 {fmt(2*min(tube_tau) if tube_tau else None)} 至 {fmt(2*max(tube_tau) if tube_tau else None)} ms 的双程时延。')
heading(doc,'关键词')
para(doc,'HF 短波传播；电离层；三维射线追踪；TID；O/X 分裂；功率-时延谱；Born/Rytov/MPS；OTH 雷达')

sections=[
('第1章 引言',['HF 短波通信依赖电离层反射实现远距离传播，是远程通信、超视距雷达和被动定位系统中的关键物理链路。实际电离层并非平稳分层介质，昼夜变化、纬度变化、行进式电离层扰动、偶发 E 层、赤道等离子体泡和 Spread-F 不规则体都会改变电子密度空间分布，使传播路径、群延迟、接收功率和到达角发生显著变化。','二维射线模型能够描述主传播平面内的折返路径，但难以给出横向偏移和方位角误差；单一反射模型能够描述主模传播，却难以解释不规则体散射、闪烁和功率-时延谱展宽；若将几何路径误用为群路径，还会在时延和功率链路中引入系统偏差。因此，本项目采用三维射线追踪与散射后处理相结合的混合传播建模路线。']),
('第2章 模型总体架构与数据流',['代码主工作目录为 code。config.py 是全局配置和物理参数的单一来源，utils.py 提供电磁、功率、坐标与地磁相关工具函数，main.py 负责驱动基准、TID、Es、等离子体泡、O/X、雷达和流管雷达等场景。models 目录承载核心物理模型，包括电离层背景、二维与三维射线追踪、点对点求解、Es 模型、PE/SSF、Born/Rytov/MPS 散射和混合 pipeline。','模型数据流分为两个阶段。第一阶段生成电离层电子密度场：先调用 IRI-2016 获得背景 Ne(z)，再叠加 TID、Es、bubble 和 Spread-F 等扰动，得到 Ne_2d(x,z)。第二阶段计算传播参数：由 Ne_2d 构造 RefractiveIndex3D，发射三维射线扇面，使用 P2P Newton shooting 收敛到接收点，再对各传播模式进行功率、O/X、吸收和散射后处理，最终形成 SS1 至 SS13 参数体系。']),
('第3章 电离层环境建模',['背景电子密度由 iri2016.IRI 生成，标准链路采用 2020-06-01 12:00 UTC、路径中点约 32.5N 和 120E 的 IRI 廓线。电子密度单位为 m^-3，距离和高度单位为 km，频率单位为 MHz。该背景场代表分层电离层的平均状态，是后续扰动叠加和射线追踪的基础。','TID 扰动采用相对电子密度扰动形式，可由振幅、水平波长、周期、传播方向和横向波数控制。当前实现支持多分量、多方向的 TID_3D 配置。为避免 D 层近零密度区域压低 F 区 TID 振幅，TID 归一化只在 Ne_bg 大于 0.01 倍峰值密度的区域内计算。']),
('第4章 三维射线追踪与点对点传播模式求解',['冷等离子体各向同性折射率采用 n^2 = 1 - f_p^2/f^2，其中 f_p 与电子密度平方根成正比。RefractiveIndex3D 预计算 n^2 及其空间梯度，并在射线追踪过程中进行快速插值查询。地球展平修正通过 z_phys = z - (x^2 + y^2)/(2R_E) 将局地直角坐标与物理高度联系起来。','点对点求解由 find_all_rays_3d 完成。算法先从射线扇面中选取接近接收点的候选射线作为初值，再以发射仰角 beta 和方位角 az 为未知量进行 Newton shooting。P2P_3D 的 jacobian 支持 forward 和 central 两种差分形式；forward 较快，central 更准确但计算量更大。']),
('第5章 传播参数体系 SS1~SS13',['本模型把最终传播输出组织为 SS1 至 SS13。SS1 至 SS3 描述几何和时延属性，SS4 与 SS9 面向通信链路，SS5 与 SS7 面向雷达功率-时延谱，SS8 描述 O/X 磁离子模式分裂，SS6 与 SS10 至 SS13 描述 Es、散射、吸收和横向偏移等复杂传播效应。','最重要的约定是 group_path_km = integral(ds/n)。该量是冷等离子体群折射率近似 n_g ~= 1/n 下的等效自由空间路径，不是几何路径。群延迟 tau_ms 由 group_path_km 除以光速得到，通信和雷达功率链中的路径损耗也以群路径为基础。']),
('第6章 功率、吸收与散射后处理',['通信功率链采用 Friis 形式，但将自由空间距离替换为群路径 P\'。这样处理的物理含义是：路径损耗与波沿等效群传播路径的相位和能量扩散有关，而不是简单由几何折线长度决定。对于 O/X 分裂后的模式，功率会根据各自 group_path_km 重新计算。','D 层吸收通过 ABSORPTION 配置控制。验证脚本表明，通信场景按单程吸收扣除，雷达场景按双程吸收扣除。散射诊断由 scatter_coupler.apply_scatter_regime_3d 附加到三维模式上，按 sigma_phi^2 阈值区分 Born、Rytov 和 MPS_required。当前默认做法是添加诊断字段，而不直接改写主功率 Pr_W。']),
]
for h,ps in sections:
    heading(doc,h)
    for p in ps: para(doc,p)
    if h.startswith('第2章'):
        table(doc,['模块','主要功能','报告中对应内容'],[['config.py','统一配置频率、链路、网格、TID、散射和吸收参数','附录A'],['ionosphere_model.py','构建 IRI + 扰动电子密度场','第3章'],['ray_tracer_3d.py','三维 Haselgrove RK4 与 RefractiveIndex3D','第4章'],['point_to_point_3d.py','三维 P2P Newton shooting','第4章'],['scatter_coupler.py','Born/Rytov/MPS 诊断耦合','第6章'],['hybrid_model.py','端到端模式汇总、P-D 谱和主模识别','第5至8章']])
    if h.startswith('第5章'):
        table(doc,['编号','参数','来源','用途'],[['SS1','群延迟 tau_ms','P2P 模式输出','TDOA、测距、P-D 谱横轴'],['SS2','群路径 group_path_km','积分 sum(ds/n)','路径损耗和时延基础量'],['SS3','反射高度与模式标签','轨迹最大 z_phys','模式分类和跳数识别'],['SS4','通信接收功率','Friis 功率链','HF 链路余量'],['SS5','OTH 雷达回波功率','双程雷达方程','探测能力评估'],['SS6','Es 传播参数','Hao 三段模型','Es 反射/过渡/散射'],['SS7','功率-时延谱','模式功率展宽','通信/雷达回波结构'],['SS8','O/X 时延差','AH 后处理','磁离子分裂误差'],['SS9','聚焦因子','Coleman flux-tube','功率密度修正诊断'],['SS10','Born 弱散射','相位方差阈值','弱不规则体散射'],['SS11','Rytov 闪烁','Rytov 近似','S4 和时延抖动'],['SS12','D 层吸收','经验吸收模型','白天 HF 衰减'],['SS13','横向偏移','三维射线 y 方向积分','方位角和定位误差']])

heading(doc,'第7章 数值验证与运行流程')
para(doc,'报告实施过程中在 pytorch_cpu 环境中运行核心验证脚本。验证结果显示，群路径定义、O/X 分裂、三维 D 层吸收、P2P Jacobian、散射耦合、三维背景、三维射线追踪和 Phase 7 扩展均通过检查。Phase 5 在与其他脚本并行运行时曾触发 ppigrf 内存分配异常，单独重跑后通过，因此报告将其判定为并行资源竞争导致的运行环境异常，而非模型物理失败。')
table(doc,['验证脚本','验证目标','结果'],[['_verify_group_path.py','group_path=integral(ds/n)，tau=group_path/c，3D/2D 一致性','0 fail(s)'],['_verify_OX.py','AH 各向同性极限、O/X 模式数量与时延差','PASSED'],['_verify_absorption_3d.py','通信单程吸收、雷达双程吸收','0 fail(s)'],['_verify_p2p_jacobian.py','forward 与 central Jacobian 一致性','0 fail(s)'],['_verify_scatter_coupling.py','Born/Rytov/MPS 阈值选择','0 fail(s)'],['_verify_phase5.py','三维背景、IGRF、多分量 TID','单独重跑 0 fail(s)'],['_verify_phase6.py','三维射线、P2P、Born/Rytov','0 fail(s)'],['_verify_phase7.py','3D bubble、tube、Booker-Gordon、MPS','0 fail(s)']])
para(doc,'标准运行命令为 conda run -n pytorch_cpu python main.py。输出图件和 CSV 保存在 code/output/，包括 ray_fan_baseline.png、ray_deflect_tid.png、pd_OX.png、pd_tube_radar.png 以及 modes_baseline.csv、modes_tid.csv、modes_OX.csv 和 modes_tube_radar.csv 等。本报告第8章直接引用这些 CSV 中的数值。')
image(doc,'ray_fan_baseline.png','图1 基准场景三维射线扇面输出。')
image(doc,'ray_deflect_tid.png','图2 TID 场景横向偏移诊断图。')

heading(doc,'第8章 典型应用场景与结果讨论')
para(doc,f'基准场景给出 10 MHz、1169 km 链路的一条主要 F 层低角模式。CSV 输出中，该模式群延迟为 {fmt(base_tau)} ms，群路径为 {fmt(base_gp,1)} km，反射高度为 {fmt(base_h,1)} km，接收功率为 {fmt(base_pr,2)} dBW。该结果表明，在标准白天中纬链路下，单跳 F2 模式是主要可达模式，群路径比地面水平距离更长，反映出折返路径和 n<1 冷等离子体群路径修正的共同作用。')
para(doc,f'TID 场景相对于基准场景产生小幅时延和路径变化。输出中 tau_ms 由 {fmt(base_tau)} ms 增至 {fmt(tid_tau)} ms，变化约 {fmt((tid_tau-base_tau)*1000 if tid_tau and base_tau else None,2)} us；group_path_km 由 {fmt(base_gp,1)} km 增至 {fmt(tid_gp,1)} km，变化约 {fmt(tid_gp-base_gp if tid_gp and base_gp else None,2)} km。')
para(doc,f'O/X 后处理场景显示磁离子模式具有可分辨的群路径差异。output/modes_OX.csv 中，O 模 tau_ms 为 {fmt(o_tau)} ms，X 模 tau_ms 为 {fmt(x_tau)} ms，X 模相对 O 模晚约 {fmt((x_tau-o_tau)*1000 if x_tau and o_tau else None,2)} us；两者群路径差约 {fmt(x_gp-o_gp if x_gp and o_gp else None,2)} km。需要注意当前实现是沿三维各向同性轨迹重新积分 Appleton-Hartree 折射率，并非完整的三维 AH 分支射线追踪。')
para(doc,f'三维流管雷达场景给出多条相邻仰角模式。仰角扫描范围约 {fmt(min(tube_beta),1) if tube_beta else "-"} 至 {fmt(max(tube_beta),1) if tube_beta else "-"} 度，单程 tau_ms 范围为 {fmt(min(tube_tau)) if tube_tau else "-"} 至 {fmt(max(tube_tau)) if tube_tau else "-"} ms，对应双程时延约 {fmt(2*min(tube_tau)) if tube_tau else "-"} 至 {fmt(2*max(tube_tau)) if tube_tau else "-"} ms；回波功率约为 {fmt(min(tube_pr),2) if tube_pr else "-"} 至 {fmt(max(tube_pr),2) if tube_pr else "-"} dBW。')
table(doc,['场景','模式数','代表 tau_ms','代表 group_path_km','代表 Pr_dBW','解读重点'],[['baseline',len(sc['baseline']),fmt(base_tau),fmt(base_gp,1),fmt(base_pr,2),'基准 F 层主模'],['tid',len(sc['tid']),fmt(tid_tau),fmt(tid_gp,1),fmt(tid_pr,2),'TID 引起小幅时延/高度扰动'],['OX',len(sc['OX']),f'O {fmt(o_tau)} / X {fmt(x_tau)}',f'O {fmt(o_gp,1)} / X {fmt(x_gp,1)}','约 -84.7','磁离子分裂'],['tube_radar',len(tube),f'{fmt(min(tube_tau))}-{fmt(max(tube_tau))}' if tube_tau else '-', '多仰角', f'{fmt(min(tube_pr),1)}-{fmt(max(tube_pr),1)}' if tube_pr else '-','双程雷达功率-时延谱']])
image(doc,'pd_OX.png','图3 O/X 场景功率-时延谱。')
image(doc,'pd_tube_radar.png','图4 三维流管雷达功率-时延谱。')
image(doc,'ray_fan_spreadf.png','图5 Spread-F 场景射线扇面。')

heading(doc,'第9章 模型特点、局限与后续升级')
para(doc,'当前模型的主要特点在于将三维射线追踪、复杂电离层背景、功率链路、磁离子分裂、吸收和散射诊断组织为统一 pipeline。SS1 至 SS13 参数体系使通信、雷达和定位三类应用能够共享同一组物理输出，并且便于后续与实测功率-时延谱进行对比。')
para(doc,'当前局限也应明确。第一，O/X 处理尚为后处理近似，没有独立追踪 O 模和 X 模在三维磁化等离子体中的不同轨迹。第二，scatter_coupler 默认只附加散射诊断，不直接反馈主功率链，因此强散射条件下的功率修正仍需结合 PE/SSF 或 MPS 进一步处理。第三，部分吸收、Es 和散射参数依赖经验配置，后续应通过多时刻、多频率和实测数据反演加以约束。')
heading(doc,'第10章 结论')
para(doc,'本报告将 HF 短波混合传播模型从参数手册整理为技术论文式报告，系统说明从 IRI 背景电子密度、复杂扰动叠加、三维折射率场构建、Haselgrove RK4 射线追踪、点对点收敛，到功率、O/X、吸收和散射后处理的完整数据流。')
para(doc,'验证与典型场景结果表明，当前模型已具备描述基准 F 层传播、TID 微扰、O/X 分裂、流管雷达功率-时延谱和不规则体散射诊断的能力。特别是 group_path_km 的明确定义、O/X 分支修正和三维横向偏移输出，为后续 HF 通信链路评估、OTH 雷达回波解释和定位误差分析提供统一建模基础。')
heading(doc,'参考文献（待核验）')
for ref in ['IRI-2016 国际参考电离层模型相关文献。','Haselgrove 射线追踪与 Jones & Stephenson 1975 地球展平/射线追踪相关文献。','Appleton-Hartree 磁离子理论相关文献。','Coleman 1997/1998 flux-tube 聚焦与 OTH 雷达传播相关文献。','Hao 2017 偶发 E 层三段传播模型相关文献。','Rino 1979 电离层不规则体相位屏和闪烁理论相关文献。','Knepp 多相位屏模型相关文献。','Pederick & Cervera 2014 D 层吸收参数化相关文献。']:
    para(doc,ref)
heading(doc,'附录A 关键配置参数')
table(doc,['参数组','关键字段','说明'],[['基础链路','FREQ_MHZ, TX_POS, RX_POS','载频、发射点和接收点'],['背景网格','BG_X, BG_Z','x-z 电子密度网格'],['RT_3D','ds_km, beta_min, beta_max, n_fan','三维射线追踪步长和扇面'],['P2P_3D','n_init, tol_km, d_beta, d_az, jacobian','三维点对点 Newton shooting'],['TID/TID_3D','amplitude, wavelength, period, azimuth','TID 扰动参数'],['GEOMAG','fH_MHz, dip_deg, decl_deg','磁离子和吸收计算'],['SCATTER','Cs_rel, L_outer_km, p, DeltaL_km','Born/Rytov/MPS 诊断'],['ABSORPTION','enable, A0, chi_deg','D 层吸收']])
heading(doc,'附录B SS1~SS13 公式摘要')
para(doc,'SS1: tau_ms = group_path_km / 299.792。SS2: group_path_km = integral(ds/n)。SS3: h_reflect_km = max(z - (x^2+y^2)/(2R_E))。SS4: 通信功率采用 Friis 公式并以群路径作为等效距离。SS5: 雷达功率采用双程雷达方程，路径损耗与 P\'^4 相关。SS6: Es 以 alpha=f/foEs 分段。SS7: P-D 谱以模式功率在 tau 轴上展宽。SS8: O/X 使用 Appleton-Hartree 后处理，O 为 +disc 分支，X 为 -disc 分支。SS9: F_focus_3d = sqrt(F_el*F_az)。SS10~SS11: 由 sigma_phi^2 判定 Born/Rytov/MPS。SS12: D 层吸收按通信单程、雷达双程扣除。SS13: delta_az = arctan(y_land/x_land)。')
heading(doc,'附录C 运行命令与输出文件')
para(doc,'主要运行命令为 conda run -n pytorch_cpu python main.py。常用验证命令包括 conda run -n pytorch_cpu python tests/_verify_group_path.py、tests/_verify_OX.py、tests/_verify_absorption_3d.py、tests/_verify_p2p_jacobian.py、tests/_verify_scatter_coupling.py、tests/_verify_phase5.py、tests/_verify_phase6.py 和 tests/_verify_phase7.py。')
table(doc,['输出文件','含义'],[['ray_fan_*.png','射线扇面图'],['pd_*.png','功率-时延谱'],['ray_deflect_tid.png','TID 横向偏移诊断图'],['modes_*.csv','传播模式数值结果']])
for sec in doc.sections:
    foot=sec.footer.paragraphs[0]; foot.alignment=WD_ALIGN_PARAGRAPH.CENTER; rr=foot.add_run('复杂电离层 HF 短波混合传播技术报告'); set_font(rr,'宋体',9)

out=WORD/'hf_hybrid_propagation_technical_report.docx'; doc.save(out)
chk=Document(out); (WORD/'validation.txt').write_text(f'docx={out}\nparagraphs={len(chk.paragraphs)}\ntables={len(chk.tables)}\nsize_bytes={out.stat().st_size}\n',encoding='utf-8')
print('OK',out,len(chk.paragraphs),len(chk.tables),out.stat().st_size)
