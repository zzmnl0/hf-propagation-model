# AH3D 射线—波动反射/散射混合模型：2021—2026 文献与实施路线讨论

- Task ID：AH3D-RAY-WAVE-ROUTE-20261104
- 日期：2026-11-04
- 项目：D:\code11\IRI01\IRI03\TID\code
- 核对版本：8fbb5e1210b0a13fb2dcfea34eee86261752bd5e
- 阶段：DISCUSSION ONLY。本文是新任务的文献和方法讨论报告，不构成实施授权；不修改原散射计划及研究代码，不启动 AH、PE、PHaRLAP、FARR 或 B1 作业。

## 1. 研究结论与项目书对应

项目书要求的主线具有可行的理论和文献基础：

1. 单程：AH3D 背景射线束 → 局部 PE／全波散射窗 → 出口模态复角谱 → AH3D 接续至 RX。
2. 双程：在同一射线—波动接口上增加 ECEF 物理流管、地面／海面散射算子和独立回程，构成 TX→电离层→地表→电离层→RX。
3. 近期已有 HF 射线—MPS 波动传播、磁化矢量 FPE、全波求解器及 OTHR 全链路回波模型。它们可以组成一条具体研究路线。此次检索尚未核实一篇同时实现“任意三维磁离子背景、射线与局部波动场双向交接、地面弥散后向散射、共同随机介质中的相干双程”的现成公开端到端模型。此处“尚未核实”不代表断言不存在。

建议采用“双层波动求解”：

- 小规模、可收敛的局部磁化全波算例先校准入射场、出口模态与功率接口。
- 在前向传播适用区域，用 PE/FPE 的 SSF／MPS 实现承担较大局部窗或多次随机实现。
- 强梯度、显著局部后向或复杂模式转换可使用局部全波／适合的双向方法。弱后向散射也可保留 DWBA 作为互补或交叉检验，但不能用它替代项目书的 PE／全波核心。

“PE/FPE、SSF、MPS”不宜当作三个互斥理论：PE/FPE 是近似方程，SSF 是分步求解策略，MPS 是相位屏和屏间波动传播框架。是否宽角、矢量以及覆盖哪些传播分量，由具体方程与传播核决定。随机叠加屏相位并不自动等于完整磁化矢量 PE。

## 2. 文献分级：单程与局部波动模块

按正式发表年份计入 2021—2026；DOI 中的年份不必等于出版年份。

### S1. Ding, M., Tong, P., Wei, Y., & Yu, L. (2021)

*Multiple Phase Screen Modeling of HF Wave Field Scintillations Caused by the Irregularities in Inhomogeneous Media*. Radio Science, 56, e2020RS007239.

- DOI：[10.1029/2020RS007239](https://doi.org/10.1029/2020RS007239)
- [出版社全文](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020RS007239)
- 分类：已有 HF 射线输入与 MPS 波动传播结合的直接参考；最适合作为单程的低成本算法入口。
- 方法相关性：使用射线和背景 Ne 作为输入，分步处理几何光学相位及屏间衍射，比较 12.7 MHz 单程天波信号的闪烁统计，并考察 TID。
- 具体边界：文中的 direct wave 经过一次 Es 电离层反射后到达接收机；不是自由空间直视信号，也不是地面散射双程。作者使用拟合介质／扰动参数、采用标量和特定谱假设，并在转向点前截断相位积分。不能把它当作任意转向区、绝对功率或完整三维 O/X 耦合的验证。

### S2. Rino, C. L., & Carrano, C. S. (2021)

*A vector theory for forward propagation in a structured ionosphere*. Journal of Atmospheric and Solar-Terrestrial Physics, 215, 105558.

- DOI：[10.1016/j.jastp.2021.105558](https://doi.org/10.1016/j.jastp.2021.105558)
- [作者公开全文](https://chuckrino.com/wp-content/uploads/2021/05/1-s2.0-S1364682621000225-main.pdf)
- 分类：磁化矢量前向波动方法，适合指导从标量 MPS/PE 升级为矢量传播。
- 可借鉴：Maxwell 方程到矢量 FPE 的构造、介电张量及分步推进。
- 边界：前向理论忽略反向贡献，推导采用均匀磁场；并非任意 ECEF 变磁场的曲线坐标耦合或真正后向散射的现成验证。向局部变磁场、双横向维和 AH 接口推广属于待验证工作。

### S3. Green, A., Longley, W. J., Oppenheim, M. M., & Young, M. A. (2025)

*An open source code for modeling radio wave propagation in earth’s ionosphere*. Frontiers in Astronomy and Space Sciences, 12, 1521497.

- DOI／[全文](https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2025.1521497/full)：10.3389/fspas.2025.1521497
- [FARR 代码](https://gitlab.com/longleywj/farr)；[文档](https://longleywj.gitlab.io/farr/)
- 分类：可复用的独立全波求解器及验证工具，不是已完成的 AH3D—散射全链路。
- 可借鉴：磁化、碰撞等离子体 FDTD，TFSF 入射与散射场组织、边界和有限带宽场输出。
- 验证口径：论文有 O/X/R/L 色散／WKB 测试；另有 PHaRLAP 对比，但该对比采用 1 MHz、二维、忽略地磁。两者不能合并声称为三维磁化 HF 射线路径的联合验证。不规则体雷达后向散射在文中属于后续应用方向。
- 本项目没有安装或运行 FARR。

### S4. Kim, E.-H., Johnson, J. R., Shiraiwa, S., Bertelli, N., Ono, M., & Wing, S. (2026)

*Effect of Density Irregularities on Radio Frequency Wave Propagation in Ionospheric Plasmas*. EPJ Web of Conferences, 346, 01008.

- DOI：[10.1051/epjconf/202634601008](https://doi.org/10.1051/epjconf/202634601008)
- [作者机构记录及摘要](https://collaborate.princeton.edu/en/publications/effect-of-density-irregularities-on-radio-frequency-wave-propagat/)
- 分类：直接相关的 HF 局部全波不规则体机制研究；2026 年发表的会议论文。
- 摘要确认：使用 Petra-M 研究 3–30 MHz 波与密度不规则体作用，包含由流体模拟得到的泡结构、多方向散射和电磁到静电模式转换。
- 可借鉴：将泡／梯度结构作为局部计算域，分析角向能量和模式转换。与参数化 bubble 路线高度相关。
- 阅读状态：本轮核对了机构摘要和出版信息；出版社全文访问失败，未核实具体网格、边界及其定量验证。不能据此断言已存在与 AH3D 的全链路耦合。

## 3. 文献分级：地表边界、双程全链路与验证

### D1. Rino, C., & Carrano, C. (2021)

*A vector theory for forward propagation in a structured ionosphere with surface reflections*. Journal of Atmospheric and Solar-Terrestrial Physics, 223, 105740.

- DOI：[10.1016/j.jastp.2021.105740](https://doi.org/10.1016/j.jastp.2021.105740)
- [作者公开稿](https://chuckrino.com/wp-content/uploads/2021/05/VectorFPESurfaceBBL.pdf)
- 分类：矢量前向传播与平滑导体表面边界的结合。
- 可借鉴：用边界积分及诱导表面源处理反射；可作为近地表波动边界方法参考。
- 边界：平滑表面反射不等于粗糙地表向雷达方向的弥散后向散射；也不自动构成去回程电离层传播模型。

### D2. Chen, Z., Ji, Y., Zhang, Y., et al. (2025)

*Skywave OTHR Full-Link Modeling and Simulation—Part I: Trans-Ionospheric Sea Clutter*. IEEE Transactions on Geoscience and Remote Sensing, 63, 5103522, 1–22.

- DOI：[10.1109/TGRS.2025.3559010](https://doi.org/10.1109/TGRS.2025.3559010)
- [作者上传的公开全文](https://www.researchgate.net/publication/390670478_Skywave_OTHR_Full-Link_Modeling_and_Simulation-Part_I_Trans-Ionospheric_Sea_Clutter)
- 分类：与双程任务最接近的近期全链路系统参考。
- 可借鉴：AH／射线计算与海面 RCS、天线、多路径、快慢时间及距离–Doppler 回波模型的整合，包含实测对照。
- 项目应用：借用其回波组织与验证观测量，以本项目的 AH＋局部 PE／全波传播算子增强其中的传播环节。
- 边界：这是射线、海面散射与回波信号的融合，不能将其称作已实现电离层不规则体 PE／全波耦合。

### D3. Chen, Z., Ji, Y., Zhang, Y., et al. (2025)

*Skywave OTHR Full-Link Modeling and Simulation—Part II: Trans-Ionospheric Multipath Target Signal*. IEEE Transactions on Radar Systems, 3, 1392–1412.

- DOI：[10.1109/TRS.2025.3623966](https://doi.org/10.1109/TRS.2025.3623966)
- [作者公开稿页面](https://www.researchgate.net/publication/396713318_Skywave_OTHR_Full-link_Modeling_and_Simulation_-_Part_II_Trans-ionospheric_Multipath_Target_Signal)
- 分类：多径目标信号与 RD 输出的辅助参考；主要基于公开摘要和作者稿检索结果。
- 可借鉴：去回多径、目标 RCS 和频散／非平稳传播对目标回波的影响。
- 边界：目标 RCS 与连续地表后向散射不是同一操作对象，不能用目标点散射替代项目书要求的地面散射。

### D4. Ruck, J. J., et al. (2024)

*On the Use of SuperDARN Ground Backscatter Measurements for Ionospheric Propagation Model Validation*. Space Weather, 22, e2024SW003916.

- DOI／[全文](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024SW003916)：10.1029/2024SW003916
- 分类：双程地面回波的观测验证路线。
- 可借鉴：地面回波前缘、群距离和仰角分布，可帮助验证射线主路径、背景 Ne 与传播模式。
- 边界：文中做过模型功率与 SNR 分布归一，不能直接当作绝对接收功率验证；也未提供项目所需的局部 PE／全波散射闭环。

补充候选：

- [Zhao 等，2026，参数化等离子体泡与斜向 HF 传播](https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2026.1753973/full)，DOI 10.3389/fspas.2026.1753973：适合参数化背景、路径和斜向电离图的测试设计；不是 PE 散射已完成的证据。
- [Nosikov 等，2025，磁各向异性介质点到点射线](https://doi.org/10.1029/2024RS008092)：可研究出口散射分量到 RX 的路径搜索方法；与散射窗本身分工不同。
- [电离层人工柱状电子云 PE 散射，2024](https://doi.org/10.1088/1402-4896/ad1959)：题名 *Electromagnetic scattering of cylindrical artificial electron cloud by using parabolic equation method*；检索摘要提及 Padé-(1,1) PE。出版社全文未获取，作者及具体角域／磁化处理未逐项核实，保留为补充候选，不作为本报告核心证据。
- [CN119557535B](https://patents.google.com/patent/CN119557535B/zh)：IRI＋TID＋宽角 PE 的技术路线候选；属于专利，不计入近年论文表。

## 4. 单程推荐路线

### 4.1 正式物理闭环

TX → AH3D 射线束 → 入射截面 → 局部矢量 PE／全波 → 出口 E/H 与模态复角谱 → AH3D 分支接续 → RX。

局部窗承担该段完整介质传播，应替换 AH 的相应传播段。若使用“背景解＋散射修正”，必须明确定义同一个窗、介质及端口下的背景算子 S0 与完整算子 S，修正为 ΔS=S−S0。背景直达部分只计一次。

按照项目书，以入射角给定方向初始条件，同时补充：

- 实际相位波矢、频率、入射时刻和局部参考坐标；群方向与相位方向分别保存。
- 有限入射束的幅度、横向波前相位／曲率和局部磁离子偏振 E/H。
- 射线束到局部截面的通量归一；有限束宽可以先作为明确输入，再由发射口径和邻近射线束约束。

出口对 E/H 进行局部传播模态分解。以切向波矢和局部 AH 色散关系求法向传播根，以外向能流筛选 O/X 分量，并保留复系数。仅在出口模态可分离、射线/WKB 有效处接回 AH；窗口边界应尽量放在不规则体之外的平滑传播区。近截止／焦散／强耦合区的交接需适合的波动表示，不能直接强分成独立射线。

接收约束不能只检查原始背景 P2P 路径：散射会形成新的 RX 可达分支，需要从出口角谱重新搜索、或建立接收端口的相应传播算子。场向／各向异性方向性保留在局部三维密度、磁场及谱中。

### 4.2 三条可行技术分支

| 分支 | 适合先回答的问题 | 对现有 AH3D 的应用 |
|---|---|---|
| 射线条件化标量 PE/SSF/MPS | 小角前向散射、闪烁、时延／角谱的低成本闭环 | 借鉴 Ding 2021；先在明示近似下验证场交接，随后增加矢量处理 |
| 曲线坐标／分扇区的矢量宽角 PE/FPE | 磁化介质中的前向衍射、偏振及较大角谱 | 借鉴 Rino–Carrano 2021；利用现有局部坐标资产，分别验证曲率项、变磁场和出口模态 |
| AH＋局部 FDTD/FEM | 强局部梯度、显著逆向散射、复杂模转换 | 借鉴 FARR 2025、Kim 2026；先做局部规范算例，再与射线端口闭环 |

这里的“局部前向”是相对于推进坐标。沿弯曲 AH 路径推进可以在地理坐标中先上升再下降；这与相对于局部推进方向的真正反向分量不同。分扇区／旋转坐标能扩大数值覆盖，但仍需证明扇区交接的幅相和通量一致。PE 名称本身不保证处理任意后向散射。

建议选择：**小全波对照确定接口 → 前向 PE 闭环并扩展计算域 → 对需要的局部结构保留全波替换。** 这既符合项目书，也不要求先实现覆盖整条千公里路径的三维全波。

## 5. 双程推荐路线

TX → 去程 AH／局部波动窗 → 地表足迹 → 地面／海面散射算子 → 回程 AH／局部波动窗 → RX。

双程与单程共享同一套局部求解器及边界端口，新增四个关键环节：

1. **ECEF 物理流管。** 由发射立体角到截面／地表面的 Jacobian 分配功率，区分垂直能流截面与地面面积。必须重新建立并验证当前 AH3D 的流管输运；legacy 的面积估计和经验聚焦因子不能直接视为已完成。
2. **地面／海面散射。** 平均功率阶段可用与面积、入／出射角及极化约定一致的 σ0。相干快慢时间和 RD 输出还需复散射幅度／散射矩阵及其时空协方差。单一 σ0 不唯一确定相位，双站也不能默认采用单站后向值。若全波域已包含地表散射，同一物理作用不再额外乘一次 σ0。
3. **独立回程。** 保留 O→地表→O、O→地表→X 等完整模式序列。总群时延为去程、地表作用及回程之和；只有经验证的对称极限才可用 2τ 和双倍吸收。磁化介质中不能未经证明使用普通互易把去程倒放；可直接计算回程，或严格建立相应伴随关系。
4. **共同介质和时间。** 去回程使用同一个 Ne(r,t)、B、碰撞及随机不规则体实现，并在各自到达时刻采样。冻结介质是可检验的时间尺度近似；不能随意使用独立随机屏，也不能把单程统计展宽机械乘二。

建议先完成“平均双程功率＋地面回波前缘／仰角分布”，再在同一物理端口上增加相干地面算子、慢时间演化和 RD 回波。两步都保留 PE／全波传播模块；区别在于观测层复杂度，而非退回经验损耗相加。

## 6. 怎样实现项目书中的损耗相加

对选定的功率归一端口或非相干传播分支，串联传输系数相乘，对应 dB 损耗相加。局部窗可定义：

T = P_selected_out / P_inc，T0 = 同窗无不规则体时的对应传输系数，

ΔL_irr = −10 log10(T/T0)。

在端口、参考面积与传播段一致时，可写为 L_total = L_background + ΔL_irr。T0=0 的新散射通道不取比值，应直接使用散射矩阵或绝对端口传输。

注意：

- 散射使能量改变方向，不全是热吸收；接收方向上的增量可正可负。
- 总 Pout/Pin≈1 并不意味着 RX 功率不变。
- 多路径先按复场或相关性规则合成，再计算可观测功率；不能把不同路径的损耗 dB 相加。
- 窗内背景吸收、几何扩散和相位只记一次。双程若已用流管与传播算子计算几何输运，就不再叠加一套群程 R^-4 或旧聚焦系数。

## 7. 真实 Ne 与参数化场景

PE／全波求解器可以共用，不必先把实测背景硬分成 bubble、spread-F 等类别。建议输入采用：

Ne_resolved(r,t) + δNe_unresolved(r,t)。

已分辨结构直接进入介电张量；未分辨结构用连续参数描述，例如扰动幅度、内外尺度、沿／横磁场相关尺度、谱斜率、漂移速度及高度包络。bubble、spread-F、Es 等作为可控测试或参数先验，而不是各自独立的求解器。

从平滑 Ne 网格不能唯一恢复这些未分辨统计量；可通过范围扫描、随机实现集合及联合观测约束，给出观测量的分布与敏感性。原有背景已包含的尺度不再重复注入随机场。

宽带复传递函数及慢时间演化可支持 P-D／功率时延相关产品、到达角、接收强度、Doppler 与模式识别。频移／Doppler 需要介质或地表的时间演化；静态介质的单频计算本身不产生真实 Doppler 展宽。散射高度、AH 转向高度与群延迟推得的虚高分别记录。

## 8. 现有代码资产与本次讨论边界

以下由 Planner／Reviewer 只读核查，绑定上述 HEAD：

- ECEF AH3D 已有路径、物理相位波矢、群时延与事件基础；没有直接提供可用于任意局部波动入射的完整复幅度、波前曲率和偏振输运。
- 当前 models/vector_pe_3d.py 有 Bishop 局部坐标、四分量复场、冷等离子体张量和通量诊断；属于研究参考。
- 既有 B1 为 blocked_by_cpu_memory_or_solver_stability／reference_uncalibrated，不能作为现成已验收 PE；本轮不重启。
- 需统一 AH／PE 的实际等离子体频率系数、时间、Ne/B、碰撞和单位。现有 PE 张量 K_FP**2 与已接受 AH 系数不同，介质采样还有 t=0.0 约定，不能默认“同介质”。
- 现有高斯初场的真空阻抗与平相位不自动适用于任意磁化入射；现有 asin(k_perp/k0) 角诊断不等于磁化出口模态映射。
- 流管模块位于 legacy；现行功率入口仍有点目标／简化双程处理。当前 P2P 的 ground 终点约束也不是任意散射截面到 RX 的现成接口。

本轮没有将上述缺口修复或加入原实施计划。

## 9. 可供后续选择的实施次序

| 次序 | 讨论中的候选工作 | 对应的最小证据 |
|---|---|---|
| A | 统一局部介质与单位功率模态端口；小规模规范全波问题 | 均匀介质 O/X、E/H 与通量；解析反射／透射；边界与网格收敛 |
| B | 单个不规则体的 AH→波动窗→AH→RX 闭环 | 零扰动恢复、移动边界不敏感、接收方向的复场／功率收敛 |
| C | PE／SSF 前向分支与局部全波对照 | 相同介质和入射场下比较角谱、幅相、时延和通量；明确后向误差 |
| D | 参数化 bubble／场向谱／spread-F 与真实 Ne 结构 | 连续参数敏感性、模式转换及实验／观测可辨识量 |
| E | ECEF 流管＋地面统计散射＋独立返回路径 | 真空雷达方程极限、面积分片收敛、去回程与模式序列；GS 前缘／仰角 |
| F | 同一时变介质的复双程回波与 RD | 地表复散射相关性、多频时延／慢时间 Doppler、自洽功率和实测对照 |

这些是路线候选及相应验证内容，不是新冻结计划、资源预算或执行批准。

最值得形成方法贡献的单元，是功率一致的 AH—波动场交接、散射出口到接收端的模态可达性，以及共同磁化介质中的相干双程。采用现有 PE／全波程序本身不宜直接宣称为首创。

## 10. 实际协作与证据

- 使用 [Codex with ChatGPT skill](C:/Users/12238/.codex/skills/codex-with-chatgpt/SKILL.md)。
- 沿用 [ChatGPT 科学讨论会话](https://chatgpt.com/c/6acefb7d-60b8-83ea-9a48-3b8c788f4da6)，新提问明确仅公开科学讨论、不恢复旧任务、不调用项目连接。已确认本轮回答完成。
- 项目代码连接尚未完成重新认证。因此 ChatGPT 的本轮作用是公开文献检索、发散分析与方法建议，没有把它表述为读取当前源码后的代码审查，也没有粘贴文件内容、diff 或日志。
- Planner 只读核对本地接口；Reviewer 独立核对关键科学表述及部分首要论文。候选路线获讨论范围 ACCEPT，最终报告是否通过以单独 review 文件为准。
- Root 独立检索和核对了上述出版信息、可获取全文或摘要，保留了访问失败。Kim、人工云论文和部分目标信号资料的阅读范围在文献条目内注明。
- 检索原始返回摘录：[web-evidence.json](../sources/ray-wave-hybrid-20261104/web-evidence.json)；摘要索引：[README.md](../sources/ray-wave-hybrid-20261104/README.md)。

本报告没有实现模型、运行数值实验或证明绝对功率预测能力；旧科学验证和 B1 门槛保持原状态。
