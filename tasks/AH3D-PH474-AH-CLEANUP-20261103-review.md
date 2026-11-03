# AH 脚本清理独立审查

**审查决定：ACCEPT，限定于本次原位逻辑归档和入口索引整理。** 本结论不表示 AH/PH 科学比对通过，也不表示旧文件已物理移动或删除。

审查对象为 `AH3D-PH474-AH-CLEANUP-20261103`，Git HEAD 为 `4be052c0c69f31a631beb8ed823c13d539036274`。工作区存在既有未提交修改；本结论同时绑定清理清单 SHA-256 `E699095F4F149A9D88FD01133827C90BB71C62098ED6DFCD07E089436D29D215` 及其记录的 13 项产物哈希。

## 阻塞项

无。

## 已核实的证据

- [清理清单](/D:/code11/IRI01/IRI03/TID/code/tasks/AH3D-PH474-AH-CLEANUP-20261103/cleanup-manifest.json)的 13 项路径唯一，文件长度及 SHA-256 均匹配。两份旧 V2 测试与[恢复快照](/D:/code11/IRI01/IRI03/TID/code/tasks/AH3D-PH474-AH-CLEANUP-20261103/rollback-manifest.json)逐字节一致，原路径保留。
- 12 项保护基线当前长度和哈希均匹配。15 项已有 tracked 改动的**路径集合**与清理前记录一致；其中未被保护基线覆盖的文件，不能仅凭路径集合断言内容逐字节未变。
- 两份 S1 清单各直接绑定 150 个根目录 Python 测试及 7 个 MATLAB 测试。独立 AST 盘点确认 38 项直接导入 `legacy.models`，其当前哈希均与两份 S1 记录一致。
- 121 项[历史脚本索引](/D:/code11/IRI01/IRI03/TID/code/tasks/AH3D-PH474-AH-CLEANUP-20261103/historical-script-index.json)所列文件存在、路径唯一，当前 SHA 均与索引一致。
- 正式轻量 runner 实际列出六项，均存在，且不包含两份旧 V2 测试。新索引对默认 `smoke`、P2P `reference` 与任务专属高精度配置的区分符合已查代码和既有证据。

本次为只读审查，未重新运行科学程序或测试；动态生成调用及仓库外调用完整性仍为 **UNKNOWN**。

## 延期事项

1. **既有 S1 版本漂移：CONFIRMED / MINOR / DEFERABLE。** 当前 `tests/_verify_p3_compatibility.py` 的 SHA 与两份 S1 记录不同；清理前该文件已为 tracked dirty，HEAD 哈希与 S1 一致。后续若重用旧 S1 原样校验，须在隔离环境恢复对应历史版本或明确处理版本差异；当前不得声称两份旧 S1 清单整体通过。
2. **部分文档链接兼容性：CONFIRMED（按普通文件路径解析）/ MINOR / DEFERABLE。** 三条含 `:行号` 的相对 Markdown 链接不能作为普通文件路径直接解析，但目标文件存在；支持行号语法的客户端是否可导航未验证。后续可在不破坏冻结清单的前提下使用兼容的链接形式。本审查文件链接在 Root 保存交接后才会存在。
3. **已有改动的内容基线：UNVERIFIED / MINOR / DEFERABLE。** 清理前 15 项 tracked 名单并非全部都有逐文件哈希；目前只能确认名单相同及 12 项受保护文件字节不变。若未来需要证明其他已有改动的精确内容未变，需查找其清理前哈希或其他独立基线。

## Planner Feedback

- **REQUIRED：** 后续引用本次结果时，将“15 项路径集合一致”与“12 项受保护文件哈希一致”分别表述；保留旧 S1 漂移和动态调用未知状态。
- **RECOMMENDED：** WGS84 主线继续使用其独立科学契约与审查，勿将此次整理或既有球面 IGRF 结果外推为一般 WGS84、一般三维或截止退化验收。
- **OPTIONAL：** 下次维护文档索引时处理行号链接兼容性；无需为此重跑科学矩阵。

---

## 已有15项内容基线的独立补核

以下为同一实际 Reviewer 在 Root 提供前包已接受基线位置后的只读补核；保留上文原始审查，不改13项冻结叶。

补核结论：可解除原审查第3条对这15项内容基线的UNVERIFIED。独立读取S5-REFLECTIGRF已接受的`identity-before.json`（02:27 UTC，早于清理基线07:29 UTC），其15条protected路径与清理前15项tracked名单完全一致；当前15/15文件的长度和SHA均匹配该基线。该身份文件及15项自身也均被281叶`science-manifest.json`精确绑定；该清单SHA为前包独审所记`52B63C220AB164E09139F365AF7444A16C0ED0E156D1A73720F221B6A37C8ED2`，身份文件SHA为`4AFF7086EF3D5F363596B2A7627E5ED0AABD40F4AABCC0B04BC43E4ADBEDB86B`。其中9项补足了本次12项基线未覆盖的tracked文件。限定ACCEPT不变；既有S1历史漂移和其他边界不变。

Reviewer最终确认：“清理前已有的15项tracked文件，当前长度和SHA-256均与前一科学包的身份基线一致，可以解除原审查中关于这15项内容基线的UNVERIFIED。该身份基线早于本次清理，并与这15项文件一起被已接受的281项科学清单绑定；其中9项补足了本次清理自身的12项保护基线未覆盖的tracked文件。原审查的限定ACCEPT不变，既有S1历史版本漂移和其他证据边界不变。”

---

记录：正文及补核来自实际 `/root/p31_sikey_review`，Root原文持久化。WGS84工程并行继续，该清理独审时间没有从WGS84执行窗口排除。
