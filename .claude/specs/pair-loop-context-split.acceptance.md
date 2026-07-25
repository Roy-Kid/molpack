---
slug: pair-loop-context-split
criteria:
  - id: ac-001
    summary: PairInputs 拥有 pair kernel 的全部读集，PackContext 不再直接持有它们
    type: code
    pass_when: |
      PackContext 上不再存在 xcart / atom_props / short_radius /
      short_radius_scale / latomnext / any_short_radius / any_fixed_atoms 字段；
      它们全部位于 PairInputs，经 `sys.inputs` 访问。
      `grep -n "pub xcart\|pub atom_props\|pub latomnext" src/context/pack_context.rs`
      的命中只出现在 PairInputs 的定义里。
    status: pending
  - id: ac-002
    summary: pair kernel 只接受读集，不再接受整个上下文
    type: code
    pass_when: |
      pair_term 与 AtomHotState::load 的签名接受 &PairInputs（或其借用形式），
      不再出现 &PackContext；src/objective.rs 中 `sys: &PackContext` 不再作为
      这两个函数的参数类型。
    status: pending
  - id: ac-003
    summary: 拆分本身不改变任何数值结果
    type: runtime
    pass_when: |
      tests/gradient.rs 全部通过（含 collective 那两条经过变异验证的）；
      tests/parallel_equivalence.rs 三条通过；
      固定 seed 下五个官方 Packmol 例子收敛且 validation 无违反。
      本条是 Task 1 的唯一功能判据——此步不动遍历，行为必须逐条不变。
    status: pending
  - id: ac-004
    summary: 五处遍历收敛为一处
    type: code
    pass_when: |
      src/objective.rs 中 `while jcart_id != NONE_IDX` 的链表遍历只出现一次
      （在 walk_chain 内）；fparc / gparc / fgparc / fparc_stats / fgparc_into
      各自只剩 sink 逻辑；牛顿第三定律的散射只写一遍。
      这是 Task 2，必须是独立于 Task 1 的 commit。
    status: pending
  - id: ac-005
    summary: 每一步都在独占节点上量过，数字进 commit message
    type: performance
    pass_when: |
      Task 1 与 Task 2 各自提交前跑过 bench_ab.sh（sbatch，独占节点），
      pair_kernel / run_iteration / objective_dispatch / restraint_eval /
      pack_end_to_end / collective_eval 对 e5ae159 基线的百分比写进对应 commit。
      门限不设——结构修正带着代价落地是允许的；但数字不得缺席。
      共享登录节点的测量不算数：它曾把 compute_f 的持平报成 +5%，
      也无法在诊断 +81% 时区分 7 µs 与 16 µs。
    status: pending
  - id: ac-006
    summary: 若 Task 2 仍然大幅回归，退回设计而不是硬扛
    type: performance
    pass_when: |
      若 Task 2 后 compute_f 相对基线劣化超过 20%（前次 PairView 尝试为 +81%），
      则不提交该形态，转而按 spec Task 3 的两条备选（sink 边界内联 /
      仅共享两个 rayon kernel）重试，并把否定结果记进 spec。
      几个百分点的代价可以带着落地；数量级的代价是设计错误的证据。
    status: pending
  - id: ac-007
    summary: 质量闸
    type: runtime
    pass_when: |
      cargo fmt --all --check、cargo clippy --all-targets --all-features -- -D warnings、
      cargo test --all-features 全部 exit 0。
      必须带 --all-features：rayon 套件是 feature-gated 的，
      普通 cargo test 会静默跳过它们——并行等价性因此曾连续四个 commit 未被验证。
    status: pending
---

# Acceptance — pair-loop-context-split

**架构优先。** 这份 spec 的交付物是结构本身，不是它带来的性能数字。
AC-005 要求每一步都测并把数字写进 commit，但**不设门限**：结构修正若代价几个
百分点，就带着这个数字落地并跟一个优化任务。AC-006 划出唯一的例外——数量级的
劣化不是代价，是设计错误的证据。

## AC-001 / AC-002 — 读集成为一个有名字的类型

缺陷是"60 字段的 bag 以共享引用穿过热路径"：这个边界没有表达 kernel 真正需要
什么，也无法与写它两个字段共存。把读集命名出来之后，累加器才留得住可变借用。

## AC-003 — 拆分不改变数值

Task 1 只搬字段，不动遍历。任何数值变化都说明搬错了东西。

## AC-004 — 遍历收敛

这是重复真正被消掉的判据。必须是独立 commit，否则回归无法二分到"拆分"还是
"共享"。

## AC-005 / AC-006 — 测量纪律

前一次尝试（291cea0）测试全绿、消掉了重复、并且 `compute_f` +81%。它是在
共享登录节点上验证的，那里读不出差别。独占节点的 A/B 是唯一算数的证据。

## AC-007 — 质量闸

`--all-features` 不是可选项，理由见条目本身。
