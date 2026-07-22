from pathlib import Path

from docx import Document
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.shared import Pt


OUT_DIR = Path(__file__).resolve().parent
DOCX_PATH = OUT_DIR / "hf_3d_reflection_scattering_model_summary.docx"
TXT_PATH = OUT_DIR / "hf_3d_reflection_scattering_model_summary.txt"


REFERENCES = [
    "Haselgrove, J. (1955). Ray theory and a new method for ray tracing. In The Physical Society (Ed.), The Physics of the Ionosphere (pp. 355-364). London: Physical Society.",
    "Haselgrove, C. B., & Haselgrove, J. (1960). Twisted ray paths in the ionosphere. Proceedings of the Physical Society, 75(3), 357-363.",
    "Jones, R. M., & Stephenson, J. J. (1975). A versatile three-dimensional ray tracing computer program for radio waves in the ionosphere. U.S. Department of Commerce, Office of Telecommunications, OT Report 75-76.",
    "Budden, K. G. (1985). The Propagation of Radio Waves: The Theory of Radio Waves of Low Power in the Ionosphere and Magnetosphere. Cambridge: Cambridge University Press.",
    "Azzarone, A., Bianchi, C., Pezzopane, M., Pietrella, M., Scotto, C., & Settimi, A. (2012). IONORT: A Windows software tool to calculate the HF ray tracing in the ionosphere. arXiv:1202.2079.",
    "Cervera, M. A., & Harris, T. J. (2014). Modeling ionospheric disturbance features in quasi-vertically incident ionograms using 3-D magnetoionic ray tracing and atmospheric gravity waves. Journal of Geophysical Research: Space Physics, 119, 431-440.",
    "Bilitza, D., Altadill, D., Truhlik, V., Shubin, V., Galkin, I., Reinisch, B., & Huang, X. (2017). International Reference Ionosphere 2016: From ionospheric climate to real-time weather predictions. Space Weather, 15(2), 418-429.",
    "Hocke, K., & Schlegel, K. (1996). A review of atmospheric gravity waves and travelling ionospheric disturbances: 1982-1995. Annales Geophysicae, 14, 917-940.",
    "Rino, C. L. (1979). A power law phase screen model for ionospheric scintillation: 1. Weak scatter. Radio Science, 14(6), 1135-1145.",
    "Rino, C. L. (1979). A power law phase screen model for ionospheric scintillation: 2. Strong scatter. Radio Science, 14(6), 1147-1155.",
    "Yeh, K. C., & Liu, C. H. (1982). Radio wave scintillations in the ionosphere. Proceedings of the IEEE, 70(4), 324-360.",
    "Knepp, D. L. (1983). Multiple phase-screen calculation of the temporal behavior of stochastic waves. Proceedings of the IEEE, 71(6), 722-737.",
    "Feit, M. D., & Fleck, J. A. (1978). Light propagation in graded-index optical fibers. Applied Optics, 17(24), 3990-3998.",
    "Barrick, D. E. (1972). First-order theory and analysis of MF/HF/VHF scatter from the sea. IEEE Transactions on Antennas and Propagation, 20(1), 2-10.",
    "Skolnik, M. I. (Ed.). (2008). Radar Handbook (3rd ed.). New York: McGraw-Hill.",
    "Coleman, C. J. (1997). On the simulation of backscatter ionograms. Journal of Atmospheric and Solar-Terrestrial Physics, 59(16), 2089-2099.",
]


SECTIONS = [
    (
        "一、目前模型总体路线",
        [
            "当前模型路线以统一三维电子密度输入为起点，目标是在给定 Ne(x,y,z,t) 的情况下计算短波传播路径、群时延、功率、P-D谱、传播模式和雷达观测量。阶段1仍使用参数化复杂电离层场作为可控输入，包括 IRI 背景、TID、Spread-F、Es 和 plasma bubble；后续接入外部三维/四维电子密度数据时，传播算法不应依赖用户预先声明扰动类型。",
            "主反射链采用三维 Appleton-Hartree Hamiltonian 射线追踪，即 AH3D。该链路在三维电子密度场和地磁场条件下计算 O/X 模式、射线轨迹、P2P传播模态、群路径和群时延。各向同性三维射线仍保留为 legacy 或性能回退路径，但研究参考输出应优先使用 AH3D，并在输出中标明 solver、wave_mode、branch_status 和 group path 定义。",
            "反射/散射耦合链路当前由 Ray-PE、散射诊断和 P-D谱核组成。Ray-PE 当前实现是 local_2p5d，即沿主射线抽取局部坐标系和二维截面，用 PE/SSF 形式给出参考诊断或可选功率修正。该实现不能等同于完整三维 PE，也不能宣称为完整三维全波散射。MPS 当前为路径相关统计核，用于相位扰动、S4、sigma_phi、sigma_tau 和 delay kernel 诊断；在完成能量闭合和文献算例验证前，不默认改变主 P-D谱。",
            "观测层区分通信单程和雷达双程。通信输出使用 tau_ms、Pr_W、Pr_dBW 和单程 P-D谱；雷达输出使用 tau_2way_ms，并区分 monostatic point target、surface clutter 和 tube radar。当前绝对功率为参考计算值，缺少实测校准、真实天线方向图、海况和接收机处理链，因此必须标记 calibration_status 或 reference_uncalibrated。",
        ],
    ),
    (
        "二、模型功能实现参考文献/理论",
        [
            "电离层背景采用 IRI/IRI-2016 作为平均电离层电子密度基础，复杂环境通过 TID、Spread-F、Es 和 bubble 参数化组件叠加。TID 的理论基础来自大气重力波和行进电离层扰动模型；Spread-F 和闪烁参考 Rino 幂律谱、Yeh-Liu 闪烁理论和多相位屏理论；Es 用薄层增强、斑块结构和 aspect sensitivity 描述；bubble 使用场向延展耗空结构描述。",
            "三维反射主链以 Appleton-Hartree 磁离子色散关系和 Haselgrove Hamiltonian 射线方程为基础。Jones-Stephenson 三维射线程序、IONORT、PHaRLAP/AU-Ray 类工具是实现和验证的重要参考。当前 AH3D first version 采用有限差分 Hamiltonian 形式，需要继续加强 near-cutoff、evanescent、分支连续性、群速度一致性和步长收敛验证。",
            "P2P传播模式求解采用发射角/方位角搜索和 Newton shooting 思路，在固定频率下寻找 TX-RX 命中的传播模式。当前实现可输出 O/X模式、反射高度、群路径、群时延和模式标签；后续需升级为扫频+扫角，以识别 MUF/LUF、模式出现/消失和跨频连续性。",
            "吸收模型采用复 Appleton-Hartree 折射率虚部路径积分思想，雷达双程按双程路径处理。当前实现已经与 AH3D 参考运行耦合，但碰撞频率、背景中性大气、电离层实测校准仍需进一步完善，legacy D-layer 吸收只能作为回退或对照。",
            "散射链路分为弱散射 Born/Rytov 判据、中强散射 MPS 和局部 PE/SSF 参考耦合。Rino、Yeh-Liu 和 Knepp 支撑相位屏和闪烁统计，PE/SSF 支撑准光学抛物方程传播，FDTD/FARR 作为未来完整三维全波散射接口与验证方向。当前代码中 local_2p5d PE 和路径相关 MPS 都是阶段性实现，不应表述为完整 3D PE 或完整 3D full-wave。",
            "雷达观测层参考 Skolnik 雷达方程、Barrick 海面 Bragg 散射和 Coleman OTH backscatter/tube tracer 思路。point target 使用 R^-4 单站雷达方程；surface clutter 和 tube radar 使用 Bragg range-Doppler projection 参考实现。由于海况、天线方向图、接收机处理链和实测校准尚未闭合，雷达功率和 P-D谱应解释为参考模型输出。",
        ],
    ),
    (
        "三、完成的工作与不足",
        [
            "当前工作已经把阶段1模型组织为“统一电子密度输入、三维反射主链、局部散射耦合、单程/双程观测输出”的框架。主要进展包括 AH3D 参考主链、复杂电离层场景分离、main2.py统一pipeline、图像背景语义修复、P-D谱语义标注和Ray-PE/MPS风险隔离。",
            "不足主要集中在四个方面。第一，ParameterizedIonosphere3DProvider 仍不能完整表达所有真实三维不规则体物理，尤其是Spread-F三维各向异性谱和Es aspect sensitivity。第二，Ray-PE仍是local_2p5d，不是完整三维PE或全波。第三，MPS是路径相关统计核，不是完整三维散射场。第四，绝对功率缺少实测校准，海况、天线方向图和接收机处理链会显著改变雷达和通信功率结果。",
        ],
    ),
    (
        "四、后续升级建议",
        [
            "下一步应优先继续收敛 main.py 到 unified_hybrid_pipeline，避免外部 Ne 接入后出现重复传播逻辑。随后应加强 AH3D 数值稳定性和O/X分支连续性，补齐 Spread-F、Es、bubble 的三维物理组件，再推进 MPS delay-angle-power kernel 能量闭合和FARR/FDTD接口验证。",
            "完整3D PE升级应作为独立阶段实施。升级路径包括定义三维传播坐标和场变量、建立AH3D射线到三维边界场的注入方式、实现双横向维度PE/SSF或FDTD局部全波域、加入PML/透明边界、做均匀介质能量守恒和文献算例对照，最后才允许影响Pr和P-D谱主结果。",
        ],
    ),
]


WORK_TABLE = [
    ("NeProvider3D与复杂电离层场", "BaseNeProvider/GridNeProvider/ParameterizedIonosphere3DProvider/ExternalNeProvider3D；main2.py通过统一provider构建baseline、TID、Spread-F、bubble、Es、full。", "统一了场景入口和图像背景；TID/Spread-F/Es/bubble背景可分离显示。", "Spread-F、Es、bubble仍为阶段性参数化；不等价于完整实测三维不规则体。"),
    ("电波单程/双程3D反射射线追踪模型：AH3D", "models/ray_tracer_ah_3d.py 与 point_to_point_ah_3d.py；main.py/main2.py参考运行启用 AH3D。", "支持O/X分支、3D轨迹、P2P模式、群路径和传播时延。", "有限差分Hamiltonian first version；near-cutoff、evanescent、步长收敛和跨频连续性仍需加强。"),
    ("P2P传播模式求解", "3D P2P固定频率扫角；AH3D模式可使用ray fan warm start。", "输出传播模式、反射高度、tau_ms、Pr和CSV语义字段。", "尚未实现扫频+扫角和跨频mode_id追踪。"),
    ("复折射率吸收", "complex_absorption_path_dB；AH3D参考运行中启用complex absorption。", "支持按路径积分吸收，雷达场景区分单程/双程。", "碰撞频率和中性大气模型仍是参考参数；缺少实测校准。"),
    ("电波传播反射/散射耦合模型：Ray-PE local_2p5d", "models/hybrid_ray_pe.py、local_frame.py、pe_propagator.py；默认不改写功率。", "建立了从射线到局部PE/SSF截面的接口和诊断字段。", "local_2p5d不是完整3D PE/全波；不能描述完整三维横向散射、路径间相干耦合和三维边界条件。"),
    ("MPS/散射P-D核", "scatter_phase_screen.py、pd_scatter_kernel.py、scatter_coupler.py；默认diagnostic/off。", "可记录sigma_phi、S4、sigma_tau、kernel_normalization和能量闭合字段。", "当前为路径相关统计核，不是完整三维散射场；scatter-to-PD默认关闭。"),
    ("P-D谱与模式输出", "hybrid_model.build_pd_spectrum、plot_pd_spectrum、save_modes_csv。", "区分通信单程、雷达双程、delay_key、pd_kernel、Ray-PE和scatter状态。", "未校准绝对功率不能解释为实测可信值。"),
    ("OTH radar、surface clutter与tube radar", "radar_clutter.py、tube_tracer_3d.py；main.py雷达场景区分point target、surface clutter、tube radar。", "支持双程时延、Bragg Doppler、range-Doppler projection参考输出。", "海况、天线方向图、接收机链、地表类型和系统损耗仍未闭合。"),
    ("main.py/main2.py运行组织", "main.py保留单场景验证；main2.py统一pipeline多场景输出到output2。", "新增RUN_PRODUCTS开关，可控制ray_fan、p2p、pd、csv等数据/图像部分。", "main.py仍有历史分支逻辑，后续应继续向unified_hybrid_pipeline收敛。"),
]


def build_txt() -> None:
    lines = ["3D复杂电离层短波反射/散射混合传播模型阶段总结", ""]
    for title, paragraphs in SECTIONS:
        lines.append(title)
        lines.extend(paragraphs)
        lines.append("")
    lines.append("模块完成情况表")
    for row in WORK_TABLE:
        lines.append("模块：{}\n当前实现：{}\n完成工作：{}\n不足/升级：{}\n".format(*row))
    lines.append("参考文献")
    lines.extend(REFERENCES)
    TXT_PATH.write_text("\n".join(lines), encoding="utf-8")


def build_docx() -> None:
    doc = Document()
    styles = doc.styles
    styles["Normal"].font.name = "SimSun"
    styles["Normal"]._element.rPr.rFonts.set(qn("w:eastAsia"), "SimSun")
    styles["Normal"].font.size = Pt(10.5)
    for style_name in ("Heading 1", "Heading 2"):
        styles[style_name].font.name = "SimHei"
        styles[style_name]._element.rPr.rFonts.set(qn("w:eastAsia"), "SimHei")

    heading = doc.add_heading("3D复杂电离层短波反射/散射混合传播模型阶段总结", 0)
    heading.alignment = WD_ALIGN_PARAGRAPH.CENTER
    subtitle = doc.add_paragraph("整理范围：当前 code 目录实现、P0-P3阶段记录、main.py/main2.py运行语义与后续升级风险。")
    subtitle.alignment = WD_ALIGN_PARAGRAPH.CENTER

    for title, paragraphs in SECTIONS:
        doc.add_heading(title, level=1)
        for text in paragraphs:
            doc.add_paragraph(text)

    doc.add_heading("模块完成情况表", level=1)
    table = doc.add_table(rows=1, cols=4)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.style = "Table Grid"
    headers = ["模块", "当前代码实现", "完成工作/输出", "不足与后续升级"]
    for idx, header in enumerate(headers):
        table.rows[0].cells[idx].text = header
    for row in WORK_TABLE:
        cells = table.add_row().cells
        for idx, text in enumerate(row):
            cells[idx].text = text
            cells[idx].vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.TOP
    for row in table.rows:
        for cell in row.cells:
            for paragraph in cell.paragraphs:
                for run in paragraph.runs:
                    run.font.size = Pt(9)

    doc.add_heading("参考文献", level=1)
    for ref in REFERENCES:
        doc.add_paragraph(ref)
    doc.save(DOCX_PATH)


if __name__ == "__main__":
    build_txt()
    build_docx()
    print(DOCX_PATH)
    print(TXT_PATH)
