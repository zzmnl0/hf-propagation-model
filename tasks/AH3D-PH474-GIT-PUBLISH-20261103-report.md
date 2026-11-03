# AH3D Git 发布范围与校验记录

Task ID：`AH3D-PH474-GIT-PUBLISH-20261103`。授权原文：“提交AH至远程git仓库”。

目标为现有仓库 `https://github.com/zzmnl0/hf-propagation-model.git`，沿用当前 `master → origin/master`。本次整理前 HEAD 为 `4be052c0c69f31a631beb8ed823c13d539036274`，已有 8 个 AH 相关本地提交尚未推送；本次发布包含这些提交及当前 AH 工作区的相关改动，采用普通快进推送，不改远程默认分支 `main`。

纳入 L2 正式生产/科研入口、SI 物理系数与吸收一致性修正、IGRF 插值精度、相关验证脚本和 AH 文档。`tasks/` 只选取总体计划及当前迁移、清理、近期反射比较的必要记录，不整体上传历史任务档案。冻结任务字节由 `.gitattributes` 保留；旧版本按既有逻辑归档方式原位保存。

MAT/NPZ 原始数组、锁文件、PH 原始反编译 C 与 native-evidence JSON、`output/`、`output2/` 生成物及未跟踪的 DOCX 草稿不随本次 AH 代码发布上传。原始数组和第三方反编译输出在本地保留；Git 检出不构成完整科学数据交付。详见[任务目录说明](./README.md)。

首次整理的 2,597 文件提交 `5d91474c85b7f4d9025733fc1b894352b44a10db` 未推送。自动审批拒绝该推送，理由是现有授权不足以确认大量任务证据可向该远程外发。该提交保留在本地分支 `codex/ah-l2-local-archive-20261103`，随后仅重组未发布的本次提交，未改写此前 8 个提交或更改冻结工作区文件；后续发布使用上述缩小范围。

本次实际检查：

- 当前 L2 候选 manifest 的 53 项文件哈希一致，漂移为 0；既有独立审查 ACCEPT 及其验证范围保持。
- 17 份新增或修改的活跃 Python 源码/测试通过 AST 解析。
- 跟踪改动通过 `git diff --check`。
- 本次纳入范围的常见凭据模式扫描无匹配；此项仅为模式检查。
- 既有 8 个未推送提交新增 blob 合计约 1.49 MiB，无超过 50 MiB 的单个 blob。

本次没有重新执行 AH/PH 科学轨迹、MATLAB、REA、完整场景或完整 P2P，也不改变物理验收、失败项和 UNKNOWN。L2 的原始回归证据见[迁移报告](./AH3D-PH474-L2-DEFAULT-20261103-report.md)和[独立审查](./AH3D-PH474-L2-DEFAULT-20261103-review.md)。本次为提交前自检，不新增独立科学审查结论。
